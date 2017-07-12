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
#define I_APP 1.81 /*uA/cm2*/ //was 1.81 for the prc stuff, experimenting with changing this, risky
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
#define TAUSYN 10			//McCarthy taui baseline = 5.0, low-dose Propofol = 10, high-dose Propofol = 20
#define USE_I_APP 0			//really should be called "USE_IAPP_STEP"
#define I_APP_START 500
#define I_APP_END 501
#define USE_LOWPROPOFOL 0	//obviously low and high propofol can't be used together, if both are 1, then lowpropofol is used
#define USE_HIGHPROPOFOL 1
#define PROPOFOL_START 1.0
#define PROPOFOL_END 100000.0
#define LOWPROP_GSYN 0.25
#define LOWPROP_TAU 10
#define HIGHPROP_GSYN 0.5
#define HIGHPROP_TAU 20
#define STARTTIME 0
#define ENDTIME 4700
#define STEPSIZE 0.01
#define DELAY 0.0			//delay must evenly divide stepsize, and it is only used if it is >= stepsize
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
#define FULLNAME "high20del0.data"
#define DBIT 1
//~ #define USE_MULTIPLIERS	0

double current[C];	//external current variable, similar to how Canavier did it
static double *del;
double gsyn, tau;
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

// Conrad code to print an array for debugging
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
		
		if (i <= 2 && USE_I_APP) {
			iapp = (time < I_APP_START || time > I_APP_END) ? I_APP : I_APP_STEP;	
		}
	
		else {
			iapp = I_APP;
		}
		
		// (((y[V + (N * i)] - (-54)) / 4) < 10e-6) ? (0.32 * 4.0) :
		current[I_NA] = G_NA * y[H + (N * i)] * pow(y[M + (N * i)], 3.0) * (y[V + (N * i)] - E_NA);
		current[I_K] =  G_K * pow(y[NV + (N * i)], 4.0) * (y[V + (N * i)] - E_K);
		current[I_M] =  G_M * y[MN + (N * i)] * (y[V + (N * i)] - E_K);
		current[I_L] =  G_L * (y[V + (N * i)] - E_L);
		
		//Ruben's crazy idea, seems to work!
		current[I_S] = 0.0;
		if (DELAY >= STEPSIZE) {
			for (j = 0; j < NN; ++j) {
				current[I_S] += weight[j + (i * NN)] * oldv[j]; //sums up products of weight and presynaptic y[S]
			}
		} else {
			for (j = 0; j < NN; ++j) {
				current[I_S] += weight[j + (i * NN)] * y[S + (N * j)]; //sums up products of weight and presynaptic y[S]
			}
		//	printf("time = %f i = %d sum = %f\n", time, i, current[I_S]);
		}		
		y[P + (N * i)] = current[I_S] ;	//sets perturbation state variable to synaptic current, doesn't affect simulation, purely for debugging purposes
		
		current[I_S] *= gsyn * (y[V + (N * i)] - E_SYN);
		
		//???????????vvvvvvv
		//~ if (USE_MULTIPLIERS) {
			//~ current[I_S] =  gsyn * ((y[S + (N * i)] * (MYCLUSTER - 1))+ (y[P + (N * i)] * (POPULATION - MYCLUSTER))) * (y[V + (N * i)] - E_SYN);
		//~ }
		//~ else {
			//~ current[I_S] = gsyn * y[S + (N * i)] * (y[V + (N * i)] - E_SYN);
		//~ }
		
		dydx[V + (N * i)] = (iapp - current[I_NA] - current[I_K] - current[I_M] - current[I_L] - current[I_S]) / CM;
		dydx[M + (N * i)] =  (( fabs(((y[V + (N * i)] + 54)) / 4) < 10e-6) ? (0.32 * 4.0) : ( f(y[V + (N * i)], 0.32, -54, 4.0) )) * (1.0 - y[M + (N * i)]) - ((fabs(((y[V + (N * i)] + 27)) / 5) < 10e-6) ? (-0.28 * -5) : ( f(y[V + (N * i)], -0.28, -27, -5.0) )) * y[M + (N * i)];
		//0.32 * (y[V + (N * i)] + 54.0) / (1.0 - exp(-(y[V + (N * i)] + 54.0) / 4.0)) * (1.0 - y[M + (N * i)]) - 0.28 * (y[V + (N * i)] + 27.0) / (exp((y[V + (N * i)] + 27.0) / 5.0) - 1.0) * y[M + (N * i)];   
		dydx[H + (N * i)] = 0.128 * exp(-(y[V + (N * i)] + 50.0) / 18.0) * (1.0 - y[H + (N * i)]) - 4.0 / (1.0 + exp(-(y[V + (N * i)] + 27.0) / 5.0)) * y[H + (N * i)];   
		dydx[NV + (N * i)] = ((fabs(((y[V + (N * i)] + 52)) / 5) < 10e-6) ? (0.032 * 5) : (f(y[V + (N * i)], 0.032, -52, 5.0) )) * (1.0 - y[NV + (N * i)]) - 0.5 * exp(-(y[V + (N * i)] + 57.0) / 40.0) * y[NV + (N * i)];
		//0.032 * (y[V + (N * i)] + 52.0) / (1.0 - exp(-(y[V + (N * i)] + 52.0) / 5.0)) * (1.0 - y[NV + (N * i)]) - 0.5 * exp(-(y[V + (N * i)] + 57.0) / 40.0) * y[NV + (N * i)];  
		dydx[MN + (N * i)] = ((fabs(((y[V + (N * i)] + 30)) / 9) < 10e-6) ? (3.209 * 0.0001 * 9.0) : (f(y[V + (N * i)], (3.209 * 0.0001), -30, 9.0) )) * (1.0 - y[MN + (N * i)]) + (((fabs(((y[V + (N * i)] + 30)) / 9) < 10e-6) ? (3.209 * 0.0001 * -9.0) : (f(y[V + (N * i)], (3.209 * 0.0001), -30, -9.0) )) * y[MN + (N * i)]);
		//above has -q in both approximation and formula, reexamine this!!!!
		//3.209 * 0.0001 * ((y[V + (N * i)] + 30.0) / (1.0 - exp(-(y[V + (N * i)] + 30.0) / 9.0)) * (1.0 - y[MN + (N * i)]) + (y[V + (N * i)] + 30.0) / (1.0 - exp((y[V + (N * i)] + 30.0) / 9.0)) * y[MN + (N * i)]); 
		
		dydx[S + (N * i)] = ((2 * (1 + tanh(y[V + (i * N)] / 4.0))) * (1 - y[S + (i * N)])) - (y[S + (i * N)] / tau);
		
/*
		synsum = 0.0;
		if (DELAY >= STEPSIZE) {
			//newest edit: *oldv should be the 2 * tanh() function of the previous voltage, moved outside of the driver to prevent waste of computational resources
			//~ for (j = 0; j < NN; ++j) {
					//~ synsum += oldv[j] * (1 - y[S + (N * i)]) * weight[i*NN+j];
			//~ }

			//?????????????????vvvvvvv
			//~ dydx[S + (N * i)] = (synsum - (y[S + (N * i)] / tau)); //uses *oldv which should be del, the delay pointer in the buffer
			
			dydx[S + (N * i)] = ((2 * (1 + tanh(y[V + (i * N)] / 4.0))) * (1 - y[S + (i * N)])) - y[S + (i * N)] / tau;
			
			//~ dydx[S + (N * i)] = 2 * (1 + tanh(*oldv / 4.0)) * (1 - y[S + (N * i)]) - y[S + (N * i)] / tau; //uses *oldv which should be del, the delay pointer in the buffer
			//~ printf("I exist at time %f.\n", time);
		}
		else {
			//~ for (j = 0; j < NN; ++j) {
					//~ synsum += (2 * (1 + tanh(y[V + (j * N)] / 4.0))) * weight[j + (i * NN)];
					//~ if (synsum >= 80.) {
						//~ fprintf(stderr, "synsum = %f time = %f y[V + (j * N)] = %f i(postsynaptic neuron) = %d j(presynaptic neuron) = %d weight[j + (i * NN)] = %f\n", synsum, time, y[V + (j * N)], i, j, weight[j + (i * NN)]);
					//~ }
			//~ }
			//~ dydx[S + (N * i)] = synsum * (NN - 1 - y[S + (N * i)])- y[S + (N * i)] / tau;
			//~ if (synsum >= 1) {
					//~ fprintf(stderr, "time = %g synsum = %g\n", time, synsum);
			//~ }
			//~ dydx[S + (N * i)] = 2 * (1 + tanh(y[V + (N * i)] / 4.0)) * (1 - y[S + (N * i)]) - y[S + (N * i)] / tau;
			
			
			
		}
*/
		
		dydx[P + (N * i)] = 0;
		
		//2*(1+tanh(y[V1 + N*j + (N * i)]/4.0))*(1-y[S1 + N*j + (N * i)])-y[S1 + N*j]/TAUSYN;  
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
	
	for (i = 0; i < n; i++) {			//first step
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
		fprintf(stderr, "Using Low Dose Propofol starting at time %f and ending at %f.\n", PROPOFOL_START, PROPOFOL_END);
	}
	if (USE_HIGHPROPOFOL) {
		fprintf(stderr, "Using High Dose Propofol starting at time %f and ending at %f.\n", PROPOFOL_START, PROPOFOL_END);
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

	derivs(time, v, dv, del, weight);
	rk4(v, dv, (N * NN), time, STEPSIZE, vout, del, weight);
	
	
	xx[0] = STARTTIME;		
	for (i = 0; i < (N * NN); i++) {
		v[i] = vout[i]; 
		y[0][i] = v[i];
		
	}

	for (k = 0; k < nstep; k++) {
		
		if (USE_LOWPROPOFOL) {	//changes gsyn to match the correct level of propofol in the simulation
			gsyn = (time > PROPOFOL_START || time < PROPOFOL_END) ? LOWPROP_GSYN: G_SYN; 
			tau = (time > PROPOFOL_START || time < PROPOFOL_END) ? LOWPROP_TAU : TAUSYN; 
		}
		else if (USE_HIGHPROPOFOL) {
			gsyn = (time > PROPOFOL_START || time < PROPOFOL_END) ? HIGHPROP_GSYN : G_SYN;
			tau = (time > PROPOFOL_START || time < PROPOFOL_END) ? HIGHPROP_TAU : TAUSYN;
		}
		else {	
			gsyn = (G_SYN);
			tau = TAUSYN;
		}
		
		del = buf[bufpos]; //moves the pointer one step ahead in the buffer
		derivs(time, v, dv, del, weight);										//actually does the step
		rk4(v, dv, (N * NN), time, STEPSIZE, vout, del, weight);		//actually does the step
		
		if (DELAY >= STEPSIZE) {
			for (j = 0; j < NN; ++j) {								//IMPORTANT: this calculates the f(vj) outside of derivs and places that into the buffer, preventing a large waste of computational time
				del[j] = vout[S + (N * j)];		//dereferences the delay pointer, and puts the previous f(vj) in the same spot
			}
		}
		if (vout[0] >= THRESHOLD && v[0] < THRESHOLD) {
			if (spikecount < (SAMPLESIZE + OFFSET)) {
				if (INTERPOLATE) {
					sptimes[spikecount] = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
					//~ vout[0] = -50.0; //fudging the interpolated value to the most recent v[0] value so that it exactly matches the spike template
					//~ fprintf(stderr, "sptimes[spikecount] with interpolation is %f. Without is %f. The difference is %f.\n", sptimes[spikecount], time, (sptimes[spikecount] - time));
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
		//~ makefullsingle(y, xx, nstep, 18, "18full.data");
		//~ makefullsingle(y, xx, nstep, 19, "19full.data");
	}
	else {
		fprintf(stderr, "\n\nSince PLONG == 0, v-n.data are not being written\n\n");
	}
	printdarr(v, N * NN);
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
