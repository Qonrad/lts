//Runge-Kutta Differential Equation Solver

//#include "nrutil.h"
//#include "nlts.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
#define I_APP 1.81 /*uA/cm2*/ 
#define I_APP_STEP 3.3
#define E_NA  50.0
#define E_K  -100.0
#define E_L  -67.0
#define E_SYN  -80.0
#define G_NA 100.0 /* mS/cm2*/
#define G_K   80.0
#define G_M   2				// Was previously almost always 2, McCarthy paper says "between 0 and 4", gmi
#define G_L   0.1
#define G_SYN  0.165		//McCarthy gi_i baseline = 0.165, low-dose Propofol = 0.25, high-dose Propofol = 0.5
#define TAUSYN 10			//McCarthy taui baseline = 5.0, low-dose Propofol = 10, high-dose Propofol = 20
#define USE_I_APP 1
#define I_APP_START 1000
#define I_APP_END 1001
#define USE_LOWPROPOFOL 1	//obviously low and high propofol can't be used together, if both are 1, then lowpropofol is used
#define USE_HIGHPROPOFOL 0
#define PROPOFOL_START 300
#define PROPOFOL_END 4000
#define LOWPROP_GSYN 0.25
#define LOWPROP_TAU 10
#define HIGHPROP_GSYN 0.5
#define HIGHPROP_TAU 20
#define STARTTIME 0
#define ENDTIME 4700
#define STEPSIZE 0.02
#define DELAY 10.0 			//delay must evenly divide stepsize, and it is only used if it is >= stepsize
#define THRESHOLD -20.0		//the voltage at which it counts a spike has occured, used to measure both nonperturbed and perturbed period for PRC
#define SAMPLESIZE 30 		//number of spikes that are averaged together to give unperturbed period
#define OFFSET 10			//number of spikes that are skipped to allow the simulation to "cool down" before it starts measuring the period
#define POPULATION 20		//number of neurons in the whole population
#define MYCLUSTER 10		//number of neurons in the simulated neuron's population

double current[C];	//external current variable, similar to how Canavier did it
static double *del;
/*
int deriv_(double *xp, double *Y, double *F) {
	extern double current[C*NN]; 
	int el,n,j;
	double time_;
	
	time_ = *xp;
	
	double iapp; //Conrad, declares iapp variable for step I_APP
	
	//printf("time = %f\n", time_);
	//printf("epsilon = %f		EPSILON = %f\n", epsilon, EPSILON);
	iapp = (time_ < 330 && time_ >300.) ? 2. : I_APP;//Conrad, controls the I_APP
	
	for(j=0;j<NN;j++)
	{
		current[I_NA1 + N*j] = G_NA*y[H]*pow(y[Mp],3.0)*(y[V1 + N*j]-E_NA);
		current[I_K1 + N*j] = G_K*pow(y[N1 + N*j],4.0)*(y[V1 + N*j] - E_K);
		current[I_M1 + N*j] = G_M*y[MN1 + N*j]*(y[V1 + N*j] - E_K);
		current[I_L1 + N*j] =  G_L*(y[V1 + N*j] - E_L);
		current[I_S1 + N*j] = G_SYN * y[S1 + N*j] * (y[V1 + N*j] - E_SYN);		//~ for (n=0;n<NN;n++) {
			//~ if(n!=j) {
				//~ current[I_S1 + N*j] =   current[I_S1 + N*j] + G_SYN*y[S1 +N*n]*(y[V1 +N*j] - E_SYN);
			//~ }
		//~ }
		F[V1 + N*j] = (iapp + EPSILON - current[I_NA1 + N*j] - current[I_K1 + N*j] - current[I_M1 + N*j] - current[I_L1 + N*j] - current[I_S1 + N*j])/CM;
		F[M1 + N*j] = 0.32*(y[V1 + N*j]+54.0)/(1.0-exp(-(y[V1 + N*j]+54.0)/4.0))*(1.0-y[M1 + N*j])-0.28*(y[V1 + N*j]+27.0)/(exp((y[V1 + N*j]+27.0)/5.0)-1.0)*y[M1 + N*j];   
		F[H1 + N*j] = 0.128*exp(-(y[V1 + N*j]+50.0)/18.0)*(1.0-y[H1 + N*j])-4.0/(1.0+exp(-(y[V1 + N*j]+27.0)/5.0))*y[H1 + N*j];   
		F[N1 + N*j] = 0.032*(y[V1 + N*j]+52.0)/(1.0-exp(-(y[V1 + N*j]+52.0)/5.0))*(1.0-y[N1 + N*j])-0.5*exp(-(y[V1 + N*j]+57.0)/40.0)*y[N1 + N*j];  
		F[MN1 + N*j] = 3.209*0.0001*((y[V1 + N*j]+30.0)/(1.0-exp(-(y[V1 + N*j]+30.0)/9.0))*(1.0-y[MN1 + N*j]) 
					  + (y[V1 + N*j]+30.0)/(1.0-exp((y[V1 + N*j]+30.0)/9.0))*y[MN1 + N*j]); 
		F[S1 + N*j] = 2*(1+tanh(y[V1 + N*j]/4.0))*(1-y[S1 + N*j])-y[S1 + N*j]/TAUSYN; //2*(1+tanh(y[V1 + N*j]/4.0))*(1-y[S1 + N*j])-y[S1 + N*j]/TAUSYN;  
	
	}
	//~ printf("current\n");
	//~ printdarr(current, C);
	return 0;
}
*/

typedef struct Templates {
	int steps;
	double *volts;
	double *init;
	double ibuf[(int)(DELAY / STEPSIZE)];
} Template;

// Conrad code to print an array for debugging
void printdarr(double *a, int numelems) {
	int i;
	printf("[");
	for (i = 0; i < numelems - 1; ++i) {
		printf("%f (%p), ", a[i], &a[i]);
	}
	printf("%f (%p)]\n", a[numelems - 1], &a[numelems - 1]);
}

inline double f(double v, double a, double th, double q) {
	return a * ((v - th) / (1 - exp(-(v - th) / q)));
}

derivs(double time, double *y, double *dydx, double *oldv) { 
	double iapp, gsyn, tau;
	extern double current[];
	
	if (USE_I_APP) {
		iapp = (time < I_APP_START || time > I_APP_END) ? I_APP : I_APP_STEP;
		
	}
	else {
		iapp = I_APP;
	}
	
	if (USE_LOWPROPOFOL) {
		gsyn = (time < PROPOFOL_START || time > PROPOFOL_END) ? (G_SYN  * (MYCLUSTER - 1)): (LOWPROP_GSYN * (MYCLUSTER - 1));
		tau = (time < PROPOFOL_START || time > PROPOFOL_END) ? TAUSYN : LOWPROP_TAU; 
	}
	else if (USE_HIGHPROPOFOL) {
		gsyn = (time < PROPOFOL_START || time > PROPOFOL_END) ? G_SYN : HIGHPROP_GSYN;
		tau = (time < PROPOFOL_START || time > PROPOFOL_END) ? TAUSYN : HIGHPROP_TAU; 
	}
	else {	
		gsyn = (G_SYN * (MYCLUSTER - 1));
		tau = TAUSYN;
	}
	
	//First cell. It is self-connected and the one that is measured for the PRC.
	//printf("%f\n", iapp);
	// (((y[V] - (-54)) / 4) < 10e-6) ? (0.32 * 4.0) :
	current[I_NA] = G_NA * y[H] * pow(y[M], 3.0) * (y[V] - E_NA);
	current[I_K] =  G_K * pow(y[NV], 4.0) * (y[V] - E_K);
	current[I_M] =  G_M * y[MN] * (y[V] - E_K);
	current[I_L] =  G_L * (y[V] - E_L);
	current[I_S] =  gsyn * (y[S] + y[P])* (y[V] - E_SYN);
	
	dydx[V] = (iapp - current[I_NA] - current[I_K] - current[I_M] - current[I_L] - current[I_S]) / CM;
	dydx[M] =  (( fabs(((y[V] + 54)) / 4) < 10e-6) ? (0.32 * 4.0) : ( f(y[V], 0.32, -54, 4.0) )) * (1.0 - y[M]) - ((fabs(((y[V] + 27)) / 5) < 10e-6) ? (-0.28 * -5) : ( f(y[V], -0.28, -27, -5.0) )) * y[M];
	//0.32 * (y[V] + 54.0) / (1.0 - exp(-(y[V] + 54.0) / 4.0)) * (1.0 - y[M]) - 0.28 * (y[V] + 27.0) / (exp((y[V] + 27.0) / 5.0) - 1.0) * y[M];   
	dydx[H] = 0.128 * exp(-(y[V] + 50.0) / 18.0) * (1.0 - y[H]) - 4.0 / (1.0 + exp(-(y[V] + 27.0) / 5.0)) * y[H];   
	dydx[NV] = ((fabs(((y[V] + 52)) / 5) < 10e-6) ? (0.032 * 5) : (f(y[V], 0.032, -52, 5.0) )) * (1.0 - y[NV]) - 0.5 * exp(-(y[V] + 57.0) / 40.0) * y[NV];
	//0.032 * (y[V] + 52.0) / (1.0 - exp(-(y[V] + 52.0) / 5.0)) * (1.0 - y[NV]) - 0.5 * exp(-(y[V] + 57.0) / 40.0) * y[NV];  
	dydx[MN] = ((fabs(((y[V] + 30)) / 9) < 10e-6) ? (3.209 * 0.0001 * 9.0) : (f(y[V], (3.209 * 0.0001), -30, 9.0) )) * (1.0 - y[MN]) + (((fabs(((y[V] + 30)) / 9) < 10e-6) ? (3.209 * 0.0001 * -9.0) : (f(y[V], (3.209 * 0.0001), -30, -9.0) )) * y[MN]);
	//above has -q in both approximation and formula, reexamine this!!!!
	//3.209 * 0.0001 * ((y[V] + 30.0) / (1.0 - exp(-(y[V] + 30.0) / 9.0)) * (1.0 - y[MN]) + (y[V] + 30.0) / (1.0 - exp((y[V] + 30.0) / 9.0)) * y[MN]); 
	if (DELAY >= STEPSIZE) {
		dydx[S] = 2 * (1 + tanh(*oldv / 4.0)) * (1 - y[S]) - y[S] / tau; //uses *oldv which should be del, the delay pointer in the buffer
	}
	else {
		dydx[S] = 2 * (1 + tanh(y[V] / 4.0)) * (1 - y[S]) - y[S] / tau;
	}
	dydx[P] = 0.0;
	//2*(1+tanh(y[V1 + N*j]/4.0))*(1-y[S1 + N*j])-y[S1 + N*j]/TAUSYN;  
}

void scan_(double *Y) {
	FILE *fopen(),*sp;
	int i;
	sp = fopen("state.data","r");
	for(i = 0; i < N; i++) {
		fscanf(sp, "%lf\n", &Y[i]);
	}
	fclose(sp);
}

void dump_(double Y[]) {
	FILE *fopen(),*sp;
	int i;
	sp = fopen("end.data","w");
	for(i=0;i<N;i++) {
		fprintf(sp,"%.16f\n",Y[i]);
	}
	fclose(sp);
}

void rk4(double y[], double dydx[], int n, double x, double h, double yout[], double *oldv) {
	int i;
	double xh, hh, h6, *dym, *dyt, *yt;
	
	dym = (double*) malloc(n * sizeof(double));
	dyt = (double*) malloc(n * sizeof(double));
	yt = (double*) malloc(n * sizeof(double));
	
	
	hh = h * 0.5;
	h6 = h / 6.0;
	xh = x + hh;
	
	for (i = 0; i < n; i++) {			//first step
		yt[i] = y[i] + hh * dydx[i];
	}
	
	derivs(xh, yt, dyt, oldv);				//second step
	
	for (i = 0; i < n; i++) {
		yt[i] = y[i] + hh * dyt[i];
	}
	
	derivs(xh, yt, dym, oldv);				//third step
	
	for (i = 0; i < n; i++) {
		yt[i] = y[i] + h * dym[i];
		dym[i] += dyt[i];
	}
	
	derivs(x + h, yt, dyt, oldv);			//fourth step
	
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


// This function is inefficient and should be optimized.
void snapshot(double** y, double *xx, int nstep, int var, double timestart, double timestop, Template *temp) {
	temp->steps = (int)((timestop - timestart) / STEPSIZE);
	temp->volts = (double*) malloc(sizeof(double) * temp->steps);
	int i, place = 0;
	for (i = 0; i < nstep + 1; ++i) {
		if (timestart < xx[i] && xx[i] < timestop) {
			temp->volts[place] = y[i][var];
			place++;
		}
	}
	return 0;
}

int main() {
	//Variables to do with simulation
	int i, k;
	double time;
	double *v, *vout, *dv;					//v = variables (state and current), vout = output variables, dv = derivatives (fstate)
	double **y, *xx; 						//results variables, y[1..N][1..NSTEP+1], xx[1..NSTEP+1]
	int nstep;								//number of steps necessary to reach ENDTIME from STARTTIME at the defined STEPSIZE
	extern double current[];				//external variable declaration
	
	//Variables to do with delay
	int dsteps = (int)(DELAY / STEPSIZE);	//number of steps in the delay (i.e. number of elements in the buffer)
	double buf[dsteps];
	del = buf;
	int bufpos;								//holds the position in the buffer that the del  pointer is at
	
	//Variables to do with initial PRC measurements
	double normalperiod;					//unperturbed period of the oscillation, difference between snd_time and fst_time;
	int psteps;								//number of steps in the unperturbed period
	double sptimes[SAMPLESIZE + OFFSET];	//array of times of spiking, differences will be averaged to find unperturbed period of oscillation
	int spikecount = 0;						//holds the location the simulation has reached in sptimes
	double spdiffs[SAMPLESIZE - 1];			//array of differences in the times of spiking, averaged to find the normalperiod
	double sumdiffs = 0;					//holds the sum of differences in times of sptimes, used for averaging
	double periodsd;						//standard deviation of the averaged periods
	
	//Variables for storing spike snapshot
	double fthresh = -1.0;					//time at which voltage crosses threshold on the way up
	double sndthresh = -1.0;				//time at which voltage crosses threshold on the way down
	Template spike;
	//~ int vsnapshotlen;						//length of voltage snapshot array
	//~ double quickshot[N];					//saves the state of state variables at the beginning of the snapshot
	//~ double bufshot[dsteps];					//saves state of the buffer at the beginning of the snapshot
	
	
	nstep = (ENDTIME - STARTTIME) / STEPSIZE;	// This assumes the entime is evenly divisible by the stepsize, which should always be true I think

	//Allocating memory for the arrays used in the calculations, maybe should switch to using stack instead of heap memory...
	y = (double**) malloc(sizeof(double*) * (nstep + 1));
	for (i = 0; i < (nstep + 1); i++) {
		y[i] = (double*) malloc(sizeof(double) * N);
	}
	xx = (double*) malloc(sizeof(double) * (nstep + 1));
	
	v = (double*) malloc(sizeof(double) * N);
	dv = (double*) malloc(sizeof(double) * N);
	vout = (double*) malloc(sizeof(double) * N);
		
	time = STARTTIME;
	scan_(v);				//scanning in initial variables (state variables only) 
	
	for (i = 0; i < (dsteps); ++i) {//sets every double in buffer to be equal to the steady state (initial) voltage that was just scanned in
		buf[i] = v[0];
	}
	
	derivs(time, v, dv, del);
	
	rk4(v, dv, N, time, STEPSIZE, vout, del);
	
	
	xx[0] = STARTTIME;		
	for (i = 0; i < N; i++) {
		v[i] = vout[i]; 
		y[0][i] = v[i];
		
	}
	for (k = 0; k < nstep; k++) {
		del = &buf[bufpos]; //moves the pointer one step ahead in the buffer
		derivs(time, v, dv, del);
		rk4(v, dv, N, time, STEPSIZE, vout, del);
		*del = vout[0];
		
		if (vout[0] >= THRESHOLD && v[0] < THRESHOLD) {
			if (spikecount < (SAMPLESIZE + OFFSET)) {
				sptimes[spikecount] = time;
			}
			++spikecount;			//incremented at the end so it can be used as position in sptimes			
		}
		
		if (spikecount > OFFSET) {	//finding time of threshold crossings for snapshot
			if (fthresh == -1.0 && vout[0] >= -50.0 && v[0] < -50.0) {
				fthresh = time;
				for (i = 0; i < dsteps; ++i) {
					spike.ibuf[i] = buf[i];
				}
			}
			else if (fthresh != -1.0 && sndthresh == -1.0 && vout[0] <= -50.0 && v[0] > -50.0) {
				sndthresh = time;
			}
		}
				
		time += STEPSIZE;
		xx[k + 1] = time;
		
		if (bufpos < dsteps - 1) {	//increments bufpos within the buffer each step
			bufpos++;
		}
		else {
			bufpos = 0;
		}
		
		for (i = 0; i < N; i++) {
			v[i] = vout[i];
			y[k + 1][i] = v[i];
		}
	}
	
	printf("This simulation counted %d spikes in all.\n", spikecount);
	if (spikecount >= (SAMPLESIZE + OFFSET)) {
		for (i = OFFSET; i < SAMPLESIZE + OFFSET - 1; ++i) {		//calculates differences between spike times to find each period
			sumdiffs += sptimes[i + 1] - sptimes[i];
			spdiffs[i - OFFSET] = sptimes[i + 1] - sptimes[i];
		}
		normalperiod = sumdiffs / SAMPLESIZE;
		psteps = normalperiod / STEPSIZE;
		periodsd = calculateSD(spdiffs, SAMPLESIZE - 1);
		//printdarr(spdiffs, SAMPLESIZE - 1);
		printf("The average unperturbed period is %f, which is approximately %d steps.\n", normalperiod, psteps);
		printf("The standard deviation is %f.\n", periodsd);
	}
	else {
		printf("There are not enough spikes to account for sample size and offset or something else has gone wrong.\n");
	}
	printf("fthresh = %f and sndthresh = %f\n", fthresh, sndthresh);
	
	snapshot(y, xx, nstep, V, fthresh, sndthresh, &spike);
	//~ vsnapshotlen = (int)((sndthresh - fthresh) / STEPSIZE);
	//~ double vsnapshot[vsnapshotlen];						//voltage snapshot from -50 to -50
	//~ snapshot(y, xx, nstep, V, fthresh, sndthresh, vsnapshot, quickshot);
	//~ printf("The snapshot contains %d steps.\n", vsnapshotlen);
	//~ printdarr(vsnapshot, vsnapshotlen);
	
	makedata(y, xx, nstep, V, "v.data");
	makedata(y, xx, nstep, M, "m.data");
	makedata(y, xx, nstep, H, "h.data");
	makedata(y, xx, nstep, NV, "n.data");
	
	//~ for (i = 0; i < (nstep + 1); i++) {		//commented out b/c it's causing errors and I don't know why :(
		//~ free(y[i]);
	//~ }
		
	//~ free(current);
	dump_(vout);
	free(spike.volts);
	free(v);
	free(dv);
	free(vout);
	//~ free(y);
	free(xx);
	
	return 0;
}
