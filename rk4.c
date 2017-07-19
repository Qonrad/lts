//Runge-Kutta Differential Equation Solver, abc

//#include "nrutil.h"
//#include "nlts.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define NN 20
#define N 7
#define C 5
#define I_NA 0  
#define I_K 1 
#define I_M 2 
#define I_L 3
#define I_S 4 
#define V 0  /* mV */
#define M 1  /* mV */
#define H 2  /* mM */
#define NV 3 
#define MN 4 
#define S 5
#define P 6
#define CM 1.00 /*uF/cm2*/
#define I_APP 2.0 /*uA/cm2*/ //was 1.81 for the prc stuff, experimenting with changing this, risky
#define I_APP_STEP -3.3
#define E_NA  50.0
#define E_K  -100.0
#define E_L  -67.0
#define E_SYN  -80.0
#define G_NA 100.0 /* mS/cm2*/
#define G_K   80.0
#define G_M   2				// Was previously almost always 2, McCarthy paper says "between 0 and 4", gmi
#define G_L   0.1
#define G_SYN  0.165		//McCarthy gi_i baseline = 0.165, low-dose Propofol = 0.25, high-dose Propofol = 0.5
#define TAUSYN 5			//McCarthy taui baseline = 5.0, low-dose Propofol = 10, high-dose Propofol = 20
#define USE_I_APP 0			//really should be called "USE_IAPP_STEP"
#define I_APP_HET 0
#define I_APP_NEURONS 0		//if I_APP enabled, directs the I_APP to affect neurons 0-I_APP_NEURONS 
#define I_APP_START 500
#define I_APP_END 501
#define USE_LOWPROPOFOL 1	
#define USE_HIGHPROPOFOL 1
#define LOW_PROPOFOL_START 0.0
#define LOW_PROPOFOL_END 2000.0
#define HIGH_PROPOFOL_START 2000.0
#define HIGH_PROPOFOL_END 5000.0
#define LOWPROP_GSYN 0.25 //should be 0.25, divided it by 20 instead of using DIVNN to exactly match Carmen's code
#define LOWPROP_TAU 10
#define HIGHPROP_GSYN 0.5
#define HIGHPROP_TAU 20
#define STARTTIME 0
#define ENDTIME 4700
#define STEPSIZE 0.01
#define DELAY 3.5			//delay must evenly divide stepsize, and it is only used if it is >= stepsize
#define THRESHOLD -50.0		//the voltage at which it counts a spike has occured, used to measure both nonperturbed and perturbed period for PRC
#define STHRESHOLD -50.0	//threshold used to measure just the spike, not the period between spikes
#define SAMPLESIZE 5 		//number of spikes that are averaged together to give unperturbed period
#define OFFSET 10			//number of spikes that are skipped to allow the simulation to "cool down" before it starts measuring the period
#define POPULATION 20		//number of neurons in the whole population
#define MYCLUSTER 10		//number of neurons in the simulated neuron's population
#define True 1
#define False 0
#define INTERPOLATE 1
#define PLONG 1
#define FULLNAME "lowhigh.data"
#define DBIT 1
#define DIVNN 1
#define G(X,Y) ( (fabs((X)/(Y))<1e-6) ? ((Y)*((X)/(Y)/2. - 1.)) : ((X)/(1. - exp( (X)/ (Y) ))) )
#define F(X,Y) ( (fabs((X)/(Y))<1e-6) ? ((Y)*(1.-(X)/(Y)/2.)) : ((X)/(exp( (X)/ (Y) ) -1)) )

double current[C];	//external current variable, similar to how Canavier did it
static double *del;
double gsyn, tau;
double iapps[NN];

void printdarr(double *a, int numelems) {
	int i;
	fprintf(stderr, "[");
	for (i = 0; i < numelems - 1; ++i) {
		fprintf(stderr, "%f (%p), ", a[i], &a[i]);
	}
	fprintf(stderr, "%f (%p)]\n", a[numelems - 1], &a[numelems - 1]);
}

static inline double f(double v, double a, double th, double q) {
	return a * ((v - th) / (1 - exp(-(v - th) / q)));
}

void derivs(double time, double *y, double *dydx, double *oldv, double* weight) { 
	double iapp, synsum;
	double *pert;
	int i, j;

	for (i = 0; i < NN; i++) {
		
		if (i >= I_APP_NEURONS && USE_I_APP) {
			iapp = (time < I_APP_START || time > I_APP_END) ? I_APP : I_APP_STEP;	
		}
	
		else {
			if (I_APP_HET) {
				iapp = iapps[i];
			}
			else {
				iapp = I_APP;
			}
		}
		
		current[I_NA] = G_NA * y[H + (N * i)] * y[M + (N * i)] * y[M + (N * i)] * y[M + (N * i)] * (y[V + (N * i)] - E_NA);
		current[I_K] =  G_K * y[NV + (N * i)] * y[NV + (N * i)] * y[NV + (N * i)] * y[NV + (N * i)] * (y[V + (N * i)] - E_K);
		current[I_M] =  G_M * y[MN + (N * i)] * (y[V + (N * i)] - E_K);
		current[I_L] =  G_L * (y[V + (N * i)] - E_L);
		
		current[I_S] = 0.0;
		if (DELAY >= STEPSIZE) {
			for (j = 0; j < NN; ++j) {
				current[I_S] += weight[j + (i * NN)] * oldv[j]; //sums up products of weight and presynaptic y[S]
			}
		} else {
			for (j = 0; j < NN; ++j) {
				current[I_S] += weight[j + (i * NN)] * y[S + (N * j)]; //sums up products of weight and presynaptic y[S]
			}
		}		
		y[P + (N * i)] = current[I_S] ;	//sets perturbation state variable to synaptic current, doesn't affect simulation, purely for debugging purposes
		
		current[I_S] *= (DIVNN) ? ((gsyn /(NN - 1)) * (y[V + (N * i)] - E_SYN)) : (gsyn * (y[V + (N * i)] - E_SYN)); //multiplies synaptic current by maximum synaptic conductance and other stuff
		
		//all of these (except for h) are using a method to prevent a divide by zero error I was encountering
		dydx[V + (N * i)] = (iapp - current[I_NA] - current[I_K] - current[I_M] - current[I_L] - current[I_S]) / CM;
		dydx[M + (N * i)] = ((0.32) * G((y[V + (N * i)] + 54), -4.) * (1. - y[M + (N * i)])) - ((0.28) * F((y[V + (N * i)] + 27), 5.) * y[M + (N * i)]);
		dydx[H + (N * i)] = 0.128 * exp(-(y[V + (N * i)] + 50.0) / 18.0) * (1.0 - y[H + (N * i)]) - 4.0 / (1.0 + exp(-(y[V + (N * i)] + 27.0) / 5.0)) * y[H + (N * i)];   
		dydx[NV + (N * i)] = ((0.032) * G((y[V + (N * i)] + 52), -5.) * (1. - y[NV + (N * i)])) - (0.5 * exp(-(y[V + (N * i)] + 57.) / 40.) * y[NV + (N * i)]); 
		dydx[MN + (N * i)] = ((3.209 * 0.0001) * G((y[V + (N * i)] + 30), -9.)  * (1.0 - y[MN + (N * i)])) + ((3.209 * 0.0001) * G((y[V + (N * i)] + 30), 9.) * y[MN + (N * i)]);
		dydx[S + (N * i)] = ((2 * (1 + tanh(y[V + (i * N)] / 4.0))) * (1 - y[S + (i * N)])) - (y[S + (i * N)] / tau);		
		dydx[P + (N * i)] = 0;
	}
	return;
}

void scan_(double *Y, int n, const char *filename) {
	FILE *fopen(),*sp;
	int i, j;
	sp = fopen(filename,"r");
	for (i = 0; i < n; ++i)
		if (fscanf(sp, "%lf\n", &Y[i]) != 1){
			fprintf(stderr, "Not enough variables in state.data file.\nNeed %d only %d given\n", N*NN, i-1);
			exit(1);
		}
	
	fclose(sp);
}


void dump_(double Y[]) {
	FILE *fopen(),*sp;
	int i, j, pos;
	sp = fopen("end.data","w");
	pos = 0;
	for (j = 0; j < NN; ++j) {
		for (i = 0; i < (N); ++i) {
			fprintf(sp, "%g\n", Y[pos]);
			++pos;
		}
	}	
	fclose(sp);
}

void rk4(double y[], double dydx[], int n, double x, double h, double yout[], double *oldv, double* weight) {
	int i;
	double xh, hh, h6, *dym, *dyt, *yt;
	
	dym = (double*) malloc(n * sizeof(double));
	dyt = (double*) malloc(n * sizeof(double));
	yt = (double*) malloc(n * sizeof(double));
	
	
	hh = h * 0.5;
	h6 = h / 6.0;
	xh = x + hh;
	
	for (i = 0; i < n; i++) {						//first step
		yt[i] = y[i] + hh * dydx[i];
	}
	
	derivs(xh, yt, dyt, oldv, weight);				//second step
	
	for (i = 0; i < n; i++) {
		yt[i] = y[i] + hh * dyt[i];
	}
	
	derivs(xh, yt, dym, oldv, weight);				//third step
	
	for (i = 0; i < n; i++) {
		yt[i] = y[i] + h * dym[i];
		dym[i] += dyt[i];
	}
	
	derivs(x + h, yt, dyt, oldv, weight);			//fourth step
	
	for (i = 0; i < n; i++) {			//Accumulate increments with proper weights.
		yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
	}
	
	free(dym);
	free(dyt);
	free(yt);
}

void makedata(double** y, double *xx, int nstep, int var, const char *filename) {	//makes a .data file with time and the specified variable, named with const char
	int i;
	FILE *fopen(),*fp;
	fp = fopen(filename, "w");
	for (i = 0; i < nstep + 1; i++) {
		fprintf(fp, "%f %f\n", xx[i], y[i][var]);
	}
	fclose(fp);
}

void makefull(double** y, double *xx, int nstep, const char *filename) {	//makes a .data file with time and the specified variable, named with const char
	int i, j;
	FILE *fopen(),*fp;
	fp = fopen(filename, "w");
	for (i = 0; i < nstep + 1; i++) {
		fprintf(fp, "%f ", xx[i]);
		for (j = 0; j < (N * NN) - 1; j++) {
			fprintf(fp, "%f ", y[i][j]);
		}
		fprintf(fp, "%f\n", y[i][(N * NN) - 1]);
	}
	fclose(fp);
}

void makeallvolts(double** y, double *xx, int nstep, const char *filename) {
	int i, j;
	FILE *fopen(),*fp;
	fp = fopen(filename, "w");
	for (i = 0; i < nstep + 1; i++) {
		fprintf(fp, "%f ", xx[i]);
		for (j = 0; j < NN - 1; j++) {
			fprintf(fp, "%f ", y[i][V + (j * N)]);
		}
		fprintf(fp, "%f\n", y[i][V + (N * (NN - 1))]);
	}
	fclose(fp);
}

void makefullsingle(double** y, double *xx, int nstep, int neuron, const char *filename) {	//makes a .data file with time and the specified variable, named with const char
	int i, j;
	FILE *fopen(),*fp;
	fp = fopen(filename, "w");
	for (i = 0; i < nstep + 1; i++) {
		fprintf(fp, "%f ", xx[i]);
		for (j = 0; j < N - 1; j++) {
			fprintf(fp, "%f ", y[i][j + (neuron * N)]);
		}
		fprintf(fp, "%f\n", y[i][(neuron * N) + (N - 1)]);
	}
	fclose(fp);
}

double calculateSD(double data[], int ndatas) {
    double sum = 0.0, mean, standardDeviation = 0.0;
	int i;
    
    for(i = 0; i < ndatas; ++i) {
        sum += data[i];
    }

    mean = sum / ndatas;

    for(i=0; i < ndatas; ++i) {
        standardDeviation += pow(data[i] - mean, 2);
	}

    return sqrt(standardDeviation / ndatas);
}

void copyab(double *a, double *b, int len) {									//function to copy array of doubles
	int i;
	for (i = 0; i < len; ++i) {
		b[i] = a[i];
	}
	return;
}
void printperiod(double *periodarray, int len, const char *filename) {	//makes a .data file with time and the specified variable, named with const char
	int i;
	FILE *fopen(),*fp;
	fp = fopen(filename, "w");
	for (i = 1; i < len + 1; i++) {
		fprintf(fp, "%d %f\n", i, periodarray[i - 1]);
	}
	fclose(fp);
}

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
	for (i = 0; i < dsteps; i++) {
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
		for (i = 0; i < (dsteps); ++i) {//sets every double in buffer(s) to be equal to the steady state (initial) voltage that was just scanned in
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
		if (USE_HIGHPROPOFOL && time > HIGH_PROPOFOL_START && time < HIGH_PROPOFOL_END) {
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
