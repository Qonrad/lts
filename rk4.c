//Runge-Kutta Differential Equation Solver, abc

//Attempting to make a unified solution to combine 2neur and fixtemplate branches
//Goal = Be able to put in a single set of parameters and run the code once to produce both a PRC and the true simulations
//Switching between branches as I have it set up currently is tedious and prone to errors. Hopefully this won't be too difficult and will make things much easier.

//Command for compiling on macbook
//gcc-7 rk4.c && ./a.out && python2 all-view2.py 20 full.data

//#include "nrutil.h"
//#include "nlts.h"
#include <stdio.h>
#include <math.h>
//#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "funcs.h"

double current[C];	//external current variable, similar to how Canavier did it
static double *del;
double gsyn, tau;
double iapps[NN];

int main() {
	gsyn = G_SYN;
	tau = TAUSYN;
	//Variables to do with simulation
	long int nstep = (ENDTIME - STARTTIME) / STEPSIZE;	// This assumes the endtime-starttime is evenly divisible by the stepsize, which should always be true I think
	fprintf(stderr, "The initial simulation will contain %d steps.\n", nstep);
	int i, j, k, pre;
	double time;
	double v[N * NN], vout[N * NN], dv[N * NN];					//v = variables (state and current), vout = output variables, dv = derivatives (fstate)
	double **y, *xx; 						//results variables, y[1..N][1..NSTEP+1], xx[1..NSTEP+1]
	extern double current[];				//external variable declaration
	double weight[NN*NN];
	scan_(weight, NN*NN, "weights.data");
	
	fprintf(stderr, "DELAY = %fms.\n", DELAY);
	if (USE_LOWPROPOFOL) {
		fprintf(stderr, "Using Low Dose Propofol starting at time %f and ending at %f.\n", LOW_PROPOFOL_START, LOW_PROPOFOL_END);
	}
	if (USE_HIGHPROPOFOL) {
		fprintf(stderr, "Using High Dose Propofol starting at time %f and ending at %f.\n", HIGH_PROPOFOL_START, HIGH_PROPOFOL_END);
	}
	fprintf(stderr, "Printing data to file ");
	fprintf(stderr, FULLNAME);
	fprintf(stderr, "\n");

	//Variables to do with delay
	int dsteps = (int)(DELAY / STEPSIZE);	//number of steps in the delay (i.e. number of elements in the buffer)
	double* buf[dsteps];
	for (i = 0; i < dsteps; ++i) {
		buf[i] = (double*) malloc(sizeof(double) * NN);		//allocates memory for an array of arrays, length depending on the number of neurons in the simulation
	}
	int bufpos = 0;								//holds the position in the buffer that the del  pointer is at
	del = buf[bufpos];
	

	//Variables to do with initial PRC measurements
	double normalperiod;					//unperturbed period of the oscillation, difference between snd_time and fst_time;
	int psteps;								//number of steps in the unperturbed period
	double sptimes[SAMPLESIZE + OFFSET];	//array of times of sp(PRCSKIP) ? ((((double)(k) - (double)(psteps)) - (double)(psteps)) / (double)(psteps)) : ((double)(k) - (double)(psteps))/ (double)(psteps);iking, differences will be averaged to find unperturbed period of oscillation
	int spikecount = 0;						//holds the location the simulation has reached in sptimes
	double spdiffs[SAMPLESIZE - 1];			//array of differences in the times of spiking, averaged to find the normalperiod
	double sumdiffs = 0;					//holds the sum of differences in times of sptimes, used for averaging
	double periodsd;						//standard deviation of the averaged periods
	
	//Allocating memory for the storage arrays, checking if I can, so that I don't run out of memory
	y = (double**) malloc(sizeof(double*) * (nstep + 1));
	for (i = 0; i < (nstep + 1); i++) {
		y[i] = (double*) malloc(sizeof(double) * (N * NN));
		if (y[i] == NULL) {
			fprintf(stderr, "Ran out of memory for storage array at y[%d]", i);
			return 0;
		}
	}
	xx = (double*) malloc(sizeof(double) * (nstep + 1));
	if (xx == NULL) {
		fprintf(stderr, "Ran out of memory for storage array xx]");
		return 0;
	} 
	
	time = STARTTIME;
	scan_(v, N*NN, "state.data");				//scanning in initial variables (state variables only) 

	if (DELAY >= STEPSIZE) {
		for (i = 0; i < (dsteps); ++i) {//sets every double in buffer(s) to be equal to 0
			for (j = 0; j < NN; ++j) {
					buf[i][j] = 0.0;
			}
		}
	}

	
	iapps[0] = I_APP;
	for (i = 1; i < NN; ++i) {						//adding heterogeneity through slightly different iapps like mccarthy does
		iapps[i] = iapps[i - 1] + 0.005;
	}
	
	derivs(time, v, dv, del, weight);
	rk4(v, dv, (N * NN), time, STEPSIZE, vout, del, weight);
	
	
	xx[0] = STARTTIME;		
	for (i = 0; i < (N * NN); ++i) {
		v[i] = vout[i]; 
		y[0][i] = v[i];
		
	}

	for (k = 0; k < nstep; ++k) {
		
		if (USE_LOWPROPOFOL && time > LOW_PROPOFOL_START && time < LOW_PROPOFOL_END) {	//changes gsyn to match the correct level of propofol in the simulation
			gsyn = (time > LOW_PROPOFOL_START && time < LOW_PROPOFOL_END) ? LOWPROP_GSYN : G_SYN; 
			tau  = (time > LOW_PROPOFOL_START && time < LOW_PROPOFOL_END) ? LOWPROP_TAU : TAUSYN; 
		}
		else if (USE_HIGHPROPOFOL && time > HIGH_PROPOFOL_START && time < HIGH_PROPOFOL_END) {
			gsyn = (time > HIGH_PROPOFOL_START && time < HIGH_PROPOFOL_END) ? HIGHPROP_GSYN : G_SYN;
			tau  = (time > HIGH_PROPOFOL_START && time < HIGH_PROPOFOL_END) ? HIGHPROP_TAU : TAUSYN;
		}
		else {	
			gsyn = (G_SYN);
			tau  = TAUSYN;
		}
		
		del = buf[bufpos]; //moves the pointer one step ahead in the buffer
		derivs(time, v, dv, del, weight);										//actually does the step
		rk4(v, dv, (N * NN), time, STEPSIZE, vout, del, weight);				//actually does the step
		
		if (DELAY >= STEPSIZE) {
			for (j = 0; j < NN; ++j) {								
				del[j] = vout[S + (N * j)];							//dereferences the delay pointer, and puts the previous y[S] in the same spot
			}
		}
		if (vout[0] >= THRESHOLD && v[0] < THRESHOLD) {
			if (spikecount < (SAMPLESIZE + OFFSET)) {
				if (INTERPOLATE) {
					sptimes[spikecount] = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
				}
				else {
					sptimes[spikecount] = time;
				}
			}
			++spikecount;			//incremented at the end so it can be used as position in sptimes			
		}
		
		time += STEPSIZE;
		xx[k + 1] = time;
		
		if (DELAY >= STEPSIZE) {
			bufpos = (++bufpos)%dsteps;
		}
		
		for (i = 0; i < (N * NN); i++) {
			v[i] = vout[i];
			y[k + 1][i] = v[i];
		}
	}
	
	if (DBIT) {
		makefull(y, xx, nstep, "full.data");
	}
	
	dump_(vout);
	if (PLONG) {
		makedata(y, xx, nstep, V, "v.data");
		makedata(y, xx, nstep, (V + N), "v2.data");
		makedata(y, xx, nstep, M, "m.data");
		makedata(y, xx, nstep, H, "h.data");
		makedata(y, xx, nstep, NV, "n.data");
		makedata(y, xx, nstep, S, "s.data");
		makeallvolts(y, xx, nstep, FULLNAME);
		
		
		//~ makefullsingle(y, xx, nstep, 0, "0full.data");
		//~ makefullsingle(y, xx, nstep, 1, "1full.data");
		//~ makefullsingle(y, xx, nstep, 2, "2full.data");
		//~ makefullsingle(y, xx, nstep, 3, "3full.data");
		//~ makefullsingle(y, xx, nstep, 4, "4full.data");
		//~ makefullsingle(y, xx, nstep, 5, "5full.data");
		//~ makefullsingle(y, xx, nstep, 6, "6full.data");
		//~ makefullsingle(y, xx, nstep, 7, "7full.data");
		//~ makefullsingle(y, xx, nstep, 8, "8full.data");
		//~ makefullsingle(y, xx, nstep, 9, "9full.data");
		//~ makefullsingle(y, xx, nstep, 10, "10full.data");
		//~ makefullsingle(y, xx, nstep, 11, "11full.data");
		//~ makefullsingle(y, xx, nstep, 12, "12full.data");
		//~ makefullsingle(y, xx, nstep, 13, "13full.data");
		//~ makefullsingle(y, xx, nstep, 14, "14full.data");
		//~ makefullsingle(y, xx, nstep, 15, "15full.data");
		//~ makefullsingle(y, xx, nstep, 16, "16full.data");
		//~ makefullsingle(y, xx, nstep, 17, "17full.data");
		makefullsingle(y, xx, nstep, 18, "18full.data");
		//~ makefullsingle(y, xx, nstep, 19, "19full.data");
	}
	else {
		fprintf(stderr, "\n\nSince PLONG == 0, v-n.data are not being written\n\n");
	}
	fprintf(stderr,"This simulation counted %d spikes of Neuron[0].\n", spikecount);
	if (spikecount >= (SAMPLESIZE + OFFSET)) {
		for (i = OFFSET; i < SAMPLESIZE + OFFSET - 1; ++i) {		//calculates differences between spike times to find each period
			sumdiffs += sptimes[i + 1] - sptimes[i];
			spdiffs[i - OFFSET] = sptimes[i + 1] - sptimes[i];
		}
		printperiod(spdiffs, SAMPLESIZE - 1, "period.data");
		normalperiod = sumdiffs / (SAMPLESIZE - 1);
		psteps = (int)round(normalperiod / STEPSIZE);
		periodsd = calculateSD(spdiffs, SAMPLESIZE - 1);
		fprintf(stderr, "The average unperturbed period is %f, which is approximately %d steps.\n", normalperiod, psteps);
		fprintf(stderr, "The standard deviation is %f.\n", periodsd);
	}
	else {
		fprintf(stderr, "There are not enough spikes to account for sample size and offset or something else has gone wrong.\n");
		fprintf(stderr, "Killing because it hasn't passed the test of spikecount >= (SAMPLESIZE + OFFSET).\n");
		fprintf(stderr, "v.data-n.data as well as end.data were still written, but no trace or prc processes occured.\n");
		return 0;
	}



	for (i = 0; i < dsteps; i++) {
		free(buf[i]);
	}
	for (i = 0; i < (nstep + 1); i++) {		
		free(y[i]);
	}
	free(xx);	
	
	return 0;
}
