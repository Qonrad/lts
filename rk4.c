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
#define USE_I_APP 1
#define I_APP_START 500
#define I_APP_END 501
#define USE_LOWPROPOFOL 0	//obviously low and high propofol can't be used together, if both are 1, then lowpropofol is used
#define USE_HIGHPROPOFOL 1
#define PROPOFOL_START 300
#define PROPOFOL_END 100000
#define LOWPROP_GSYN 0.25
#define LOWPROP_TAU 10
#define HIGHPROP_GSYN 0.5
#define HIGHPROP_TAU 20
#define STARTTIME 0
#define ENDTIME 4700
#define STEPSIZE 0.01
#define DELAY 0.0 			//delay must evenly divide stepsize, and it is only used if it is >= stepsize
#define THRESHOLD -50.0		//the voltage at which it counts a spike has occured, used to measure both nonperturbed and perturbed period for PRC
#define STHRESHOLD -50.0	//threshold used to measure just the spike, not the period between spikes
#define SAMPLESIZE 5 		//number of spikes that are averaged together to give unperturbed period
#define OFFSET 20			//number of spikes that are skipped to allow the simulation to "cool down" before it starts measuring the period
#define POPULATION 20		//number of neurons in the whole population
#define MYCLUSTER 10		//number of neurons in the simulated neuron's population
#define DO_PRC 1			//toggle for prc
#define DO_TRACE 1			//toggles doing trace for a single (or multiple phase perturbations) but each is recorded individually
#define TPHASE 0.985
#define INTERVAL 100			//number of intervals prc analysis will be done on
#define True 1
#define False 0
#define INTERPOLATE 1
#define PLONG 1

double current[C];	//external current variable, similar to how Canavier did it
static double *del;
static int prcmode;
static int pertmode;
double *pert;
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

typedef struct Phipairs {
	double phase;
	double fphi1;
	double fphi2;
} Phipair;

typedef struct Templates {
	int steps;
	double *volts;
	double init[N];
	double ibuf[(int)(DELAY / STEPSIZE)];
	int bufpos;
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

void derivs(double time, double *y, double *dydx, double *oldv) { 
	double iapp, gsyn, tau;
	extern double *pert;
	
	
	if (USE_I_APP && !(prcmode)) {
		iapp = (time < I_APP_START || time > I_APP_END) ? I_APP : I_APP_STEP;
		
	}
	else {
		iapp = I_APP;
	}
	
	if (USE_LOWPROPOFOL) {
		gsyn = (time < PROPOFOL_START || time > PROPOFOL_END || prcmode) ? (G_SYN): (LOWPROP_GSYN);
		tau = (time < PROPOFOL_START || time > PROPOFOL_END || prcmode) ? TAUSYN : LOWPROP_TAU; 
	}
	else if (USE_HIGHPROPOFOL) {
		gsyn = (time < PROPOFOL_START || time > PROPOFOL_END || prcmode) ? G_SYN : HIGHPROP_GSYN;
		tau = (time < PROPOFOL_START || time > PROPOFOL_END || prcmode) ? TAUSYN : HIGHPROP_TAU;
		//~ printf("%d", prcmode);
		//~ if (prcmode && pertmode) {
			//~ fprintf(stderr, "pertmode = %d, time = %f\n", pertmode, time);
		//~ } 
	}
	else {	
		gsyn = (G_SYN);
		tau = TAUSYN;
	}
	
	//cell is self-connected and the one that is measured for the PRC.
	// (((y[V] - (-54)) / 4) < 10e-6) ? (0.32 * 4.0) :
	current[I_NA] = G_NA * y[H] * pow(y[M], 3.0) * (y[V] - E_NA);
	current[I_K] =  G_K * pow(y[NV], 4.0) * (y[V] - E_K);
	current[I_M] =  G_M * y[MN] * (y[V] - E_K);
	current[I_L] =  G_L * (y[V] - E_L);
	current[I_S] =  gsyn * ((y[S] * (MYCLUSTER - 1))+ (y[P] * (POPULATION - MYCLUSTER))) * (y[V] - E_SYN);
	
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
	if (pertmode) {
		dydx[P] = 2 * (1 + tanh(*pert / 4.0)) * (1 - y[P]) - y[P] / tau;	//should probably be Ps instead of Ss
	}
	else {
		dydx[P] = 2 * (1 + tanh((THRESHOLD) / 4.0)) * (1 - y[P]) - y[P] / tau;
	}
	//2*(1+tanh(y[V1 + N*j]/4.0))*(1-y[S1 + N*j])-y[S1 + N*j]/TAUSYN;  
	return;
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

void makefull(double** y, double *xx, int nstep, const char *filename) {	//makes a .data file with time and the specified variable, named with const char
	int i, j;
	FILE *fopen(),*fp;
	fp = fopen(filename, "w");
	for (i = 0; i < nstep + 1; i++) {
		fprintf(fp, "%f ", xx[i]);
		for (j = 0; j < N - 1; j++) {
			fprintf(fp, "%f ", y[i][j]);
		}
		fprintf(fp, "%f\n", y[i][N - 1]);
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


// This function is inefficient and should be optimized.
void snapshot(double** y, double *xx, int nstep, int var, double timestart, double timestop, Template *temp) {
	temp->steps = (int)(round((timestop - timestart) / STEPSIZE));
	temp->volts = (double*) malloc(sizeof(double) * temp->steps);
	int i, place = 0;
	if (var == V) {
		for (i = 1; i < nstep + 1; ++i) {
			if (timestart < xx[i] && xx[i] < timestop) {
				temp->volts[place] = y[i][var];
				place++;
			}
		}
		temp->volts[0] = STHRESHOLD;
		temp->volts[temp->steps - 1] = STHRESHOLD;
	}
	else {
		for (i = 0; i < nstep + 1; ++i) {
		if (timestart < xx[i] && xx[i] < timestop) {
			temp->volts[place] = y[i][var];
			place++;
		}
	}
	}
	return;
}

//Prints all the information about a completed template
void printemp(Template *temp) {//will cause an error if the template doesn't have everything in it
	printf("This template contains %d steps.\n", temp->steps);
	printf("Array of Voltages\n");
	printdarr(temp->volts, temp->steps);
	printf("Initial array of state variables\n");
	printdarr(temp->init, N);
	printf("Initial array of buffer\n");
	printdarr(temp->ibuf, (int)(DELAY / STEPSIZE));
	
}

void printphi(Phipair *p, int interval, int f, const char *filename) {	//makes a .data file with the prc
	int i;
	FILE *fopen(),*fp;
	fp = fopen(filename, "w");
	if (f == 1) {
		for (i = 0; i < interval + 1; i++) {
			fprintf(fp, "%f %f\n", p[i].phase, p[i].fphi1);
		}
	}
	else if (f == 2) {
		for (i = 0; i < interval + 1; i++) {
			fprintf(fp, "%f %f\n", p[i].phase, p[i].fphi2);
		}
	}
	fclose(fp);
}
void unpertsim(double normalperiod, Template spike, int tracedata, const char* tracename) {
	int i, k, nstep, psteps;
	double time;
	double v[N], vout[N], dv[N];			//v = variables (state and current), vout = output variables, dv = derivatives (fstate)
	double **y, *xx; 						//results variables, y[1..N][1..NSTEP+1], xx[1..NSTEP+1]
	extern double current[];				//external variable declaration
	
	//Variables to do with delay
	int dsteps = (int)(DELAY / STEPSIZE);	//number of steps in the delay (i.e. number of elements in the buffer)
	double buf[dsteps];
	del = buf;
	int bufpos;
	
	nstep = (int)round((5.0 * normalperiod) / STEPSIZE);
	psteps = (int)round(normalperiod / STEPSIZE);

	//Allocating memory for the storage arrays, checking if I can, so that I don't run out of memory
	y = (double**) malloc(sizeof(double*) * (nstep + 1));
	for (i = 0; i < (nstep + 1); i++) {
		y[i] = (double*) malloc(sizeof(double) * N);
		if (y[i] == NULL) {
			fprintf(stderr, "Ran out of memory for storage array at y[%d]", i);
			return;
		}
	}
	xx = (double*) malloc(sizeof(double) * (nstep + 1));
	if (xx == NULL) {
		fprintf(stderr, "Ran out of memory for storage array xx]");
		return;
	}
	
	//setting variables to necessary initial values
	time = 0;
	prcmode = True;
	for (i = 0; i < N; ++i) {
			dv[i] = 0.0;
			v[i] = 0.0;
			vout[i] = 0.0;
			current[i] = 0.0;
	}
	
	copyab(spike.init, v, N);
	copyab(spike.ibuf, buf, (int)(DELAY / STEPSIZE));
	bufpos = spike.bufpos;
	del = &buf[bufpos]; //moves the pointer to the correct initial position in the buffer, unnecessary i believe
	
	derivs(time, v, dv, del);	//does running this one effect the phase/ perturbation? I don't think so but I'm not sure.
	rk4(v, dv, N, time, STEPSIZE, vout, del);
	
	xx[0] = time;
	for (i = 0; i < N; i++) {
		v[i] = vout[i]; 
		y[0][i] = v[i];
	}
	
	for (k = 0; k < nstep; k++) {
		del = &buf[bufpos]; //moves the pointer one step ahead in the buffer
		derivs(time, v, dv, del);
		rk4(v, dv, N, time, STEPSIZE, vout, del);
		*del = vout[0];
				
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
	printf("The trace was unperturbed. psteps = %d (%f ms)\n", psteps, (double)psteps * STEPSIZE);
	
	
	if (tracedata) {
		char a[(int)strlen(tracename) + (int)strlen("*___.data")];
		sprintf(a, "%sv.data", tracename);
		makedata(y, xx, nstep, V, a);
		sprintf(a, "%sm.data", tracename);	
		makedata(y, xx, nstep, M, a);
		sprintf(a, "%sh.data", tracename);
		makedata(y, xx, nstep, H, a);
		sprintf(a, "%sn.data", tracename);
		makedata(y, xx, nstep, NV, a);
		sprintf(a, "%sl.data", tracename);
		makedata(y, xx, nstep, MN, a);
		sprintf(a, "%ss.data", tracename);
		makedata(y, xx, nstep, S, a);
		sprintf(a, "%sp.data", tracename);
		makedata(y, xx, nstep, P, a);
		sprintf(a, "%sfull.data", tracename);
		makefull(y, xx, nstep, a);
	}
	for (i = 0; i < (nstep + 1); i++) {		
		free(y[i]);
	}
	free(xx);
	
}

void pertsim(double normalperiod, Template spike, Phipair *trace, int tracedata, const char* tracename) {
	int i, k, nstep, targstep, flag, pertpos, psteps;
	double time, flag1, flag2;
	double v[N], vout[N], dv[N];			//v = variables (state and current), vout = output variables, dv = derivatives (fstate)
	double **y, *xx; 						//results variables, y[1..N][1..NSTEP+1], xx[1..NSTEP+1]
	extern double current[];				//external variable declaration
	
	//Variables to do with delay
	int dsteps = (int)(DELAY / STEPSIZE);	//number of steps in the delay (i.e. number of elements in the buffer)
	double buf[dsteps];
	del = buf;
	int bufpos;
	
	psteps = (int)round(normalperiod / STEPSIZE);
	nstep = (int)round((5.0 * normalperiod) / STEPSIZE);
	
	//Allocating memory for the storage arrays, checking if I can, so that I don't run out of memory
	y = (double**) malloc(sizeof(double*) * (nstep + 1));
	for (i = 0; i < (nstep + 1); i++) {
		y[i] = (double*) malloc(sizeof(double) * N);
		if (y[i] == NULL) {
			fprintf(stderr, "Ran out of memory for storage array at y[%d]", i);
			return;
		}
	}
	xx = (double*) malloc(sizeof(double) * (nstep + 1));
	if (xx == NULL) {
		fprintf(stderr, "Ran out of memory for storage array xx]");
		return;
	}
	
	//setting variables to necessary initial values
	time = 0;
	prcmode = True;
	targstep = (int)round((double)psteps * trace->phase);
	flag1 = 0.0;
	flag2 = 0.0;
	for (i = 0; i < N; ++i) {
			dv[i] = 0.0;
			v[i] = 0.0;
			vout[i] = 0.0;
			current[i] = 0.0;
	}
	
	//~ printf("\n\n\n\nAttempting to use template within pertsim!\n");
	//~ printemp(&spike);
	copyab(spike.init, v, N);
	copyab(spike.ibuf, buf, (int)(DELAY / STEPSIZE));
	bufpos = spike.bufpos;
	del = &buf[bufpos]; //moves the pointer to the correct initial position in the buffer, unnecessary i believe
	pert = spike.volts;
	pertpos = 0;
	//~ printf("Template successfully initiated!\n");
	
	derivs(time, v, dv, del);	//does running this one effect the phase/ perturbation? I don't think so but I'm not sure.
	rk4(v, dv, N, time, STEPSIZE, vout, del);
	
	xx[0] = time;
	for (i = 0; i < N; i++) {
		v[i] = vout[i]; 
		y[0][i] = v[i];
	}
	
	//~ printf("\nBeginning main for loop!\n\n");
	for (k = 0; k < nstep; k++) {			
			if (k == targstep && trace->phase >= 0.0) {	//activates perturbation mode if on correct step, allows derivs() to start using the "perturbation synapse" (a [pre-recorded stimulus of the same identical neuron)
				printf("k = %d\n", k);
				pertmode = True;
				flag = True;
			}
			del = &buf[bufpos]; //moves the pointer one step ahead in the buffer
			derivs(time, v, dv, del);
			rk4(v, dv, N, time, STEPSIZE, vout, del);
			*del = vout[0];
				
			time += STEPSIZE;
			xx[k + 1] = time;
		
			if (bufpos < dsteps - 1) {	//increments bufpos within the buffer each step
				bufpos++;
			}
			else {
				bufpos = 0;
			}
			if (pertmode) {
				if (pertpos < spike.steps - 1) {
					pertpos++;
					pert = &spike.volts[pertpos];
				}
				else {
					pertmode = False;
				}
			}
			if (flag && vout[0] >= THRESHOLD && v[0] < THRESHOLD && flag1 == 0.0 && trace->phase >= 0.0) {
				
				if (INTERPOLATE) {
					flag1 = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
					trace->fphi1 = ((double)(flag1) - (double)(normalperiod)) / (double)(normalperiod);
					printf("The trace was done at phase %f. The step targeted was %d (time of flag1= %f), and the total number of steps in the unperturbed period was %d (%f ms).\n", trace->phase, targstep, time, psteps, normalperiod);
					printf("f(phi)1 is %f.\n", trace->fphi1);
				}
				else {
					flag1 = k;
					trace->fphi1 = ((double)(k) - (double)(psteps))/ (double)(psteps);
				}
			}
			else if (flag && vout[0] >= THRESHOLD && v[0] < THRESHOLD && flag1 != 0.0 && flag2 == 0.0) {
				if (INTERPOLATE) {
					flag2 = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
					trace->fphi2 = ((double)(flag2) - (double)(flag1) - (double)(normalperiod))/ (double)(normalperiod);
				}
				else {
					flag2 = k;
					trace->fphi2 = ((double)(k) - (double)(flag1))/ (double)(psteps);
					//~ printf("k %d ptime %f targstep %d prc[prcpos].phase %f time %f prc[prcpos].fphi1 %f flag1 %d flag2 %d prc[prcpos].fphi2 %f\n", k, ptime, targstep, prc[prcpos].phase, time, prc[prcpos].fphi1, flag1, flag2, prc[prcpos].fphi2);
				}
			}
			for (i = 0; i < N; i++) {
				v[i] = vout[i];
				y[k + 1][i] = v[i];
			}
		}
	printf("targstep = %d\n", targstep);
	if (tracedata) {
		char a[(int)strlen(tracename) + (int)strlen("*___.data")];
		sprintf(a, "%sv.data", tracename);
		makedata(y, xx, nstep, V, a);
		sprintf(a, "%sm.data", tracename);	
		makedata(y, xx, nstep, M, a);
		sprintf(a, "%sh.data", tracename);
		makedata(y, xx, nstep, H, a);
		sprintf(a, "%sn.data", tracename);
		makedata(y, xx, nstep, NV, a);
		sprintf(a, "%sl.data", tracename);
		makedata(y, xx, nstep, MN, a);
		sprintf(a, "%ss.data", tracename);
		makedata(y, xx, nstep, S, a);
		sprintf(a, "%sp.data", tracename);
		makedata(y, xx, nstep, P, a);
		sprintf(a, "%sfull.data", tracename);
		makefull(y, xx, nstep, a);
	}
	for (i = 0; i < (nstep + 1); i++) {		
		free(y[i]);
	}
	free(xx);	
}

void prc(Template spike, int interval, double normalperiod) {
	Phipair prc [interval + 1];
	int i;
	for (i = 0; i < interval + 1; ++i) {
		prc[i].phase = ((double)i) * (1.0 / (double)interval);
		pertsim(normalperiod, spike, &prc[i], 0, "test");
	}
	printphi(prc, interval, 1, "prc1.data");
	printphi(prc, interval, 2, "prc2.data");
}

int main() {
	//Variables to do with simulation
	long int nstep = (ENDTIME - STARTTIME) / STEPSIZE;	// This assumes the endtime-starttime is evenly divisible by the stepsize, which should always be true I think
	printf("The initial simulation will contain %d steps.\n", nstep);
	int i, k;
	double time;
	double v[N], vout[N], dv[N];					//v = variables (state and current), vout = output variables, dv = derivatives (fstate)
	double **y, *xx; 						//results variables, y[1..N][1..NSTEP+1], xx[1..NSTEP+1]
	extern double current[];				//external variable declaration
	
	//Variables to do with delay
	int dsteps = (int)(DELAY / STEPSIZE);	//number of steps in the delay (i.e. number of elements in the buffer)
	double buf[dsteps];
	del = buf;
	int bufpos;								//holds the position in the buffer that the del  pointer is at
	
	//Variables to do with initial PRC measurements
	double normalperiod;					//unperturbed period of the oscillation, difference between snd_time and fst_time;
	int psteps;								//number of steps in the unperturbed period
	double sptimes[SAMPLESIZE + OFFSET];	//array of times of sp(PRCSKIP) ? ((((double)(k) - (double)(psteps)) - (double)(psteps)) / (double)(psteps)) : ((double)(k) - (double)(psteps))/ (double)(psteps);iking, differences will be averaged to find unperturbed period of oscillation
	int spikecount = 0;						//holds the location the simulation has reached in sptimes
	double spdiffs[SAMPLESIZE - 1];			//array of differences in the times of spiking, averaged to find the normalperiod
	double sumdiffs = 0;					//holds the sum of differences in times of sptimes, used for averaging
	double periodsd;						//standard deviation of the averaged periods
	
	//Variables for storing spike snapshot
	double fthresh = -1.0;					//time at which voltage crosses threshold on the way up
	double sndthresh = -1.0;				//time at which voltage crosses threshold on the way down
	Template spike;

	
	//Allocating memory for the storage arrays, checking if I can, so that I don't run out of memory
	y = (double**) malloc(sizeof(double*) * (nstep + 1));
	for (i = 0; i < (nstep + 1); i++) {
		y[i] = (double*) malloc(sizeof(double) * N);
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
	
	//Variables to do with conducting prc
	extern int pertmode;
	extern double *pert;
	
	
	
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
				if (INTERPOLATE) {
					sptimes[spikecount] = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
					//~ printf("sptimes[spikecount] with interpolation is %f. Without is %f. The difference is %f.\n", sptimes[spikecount], time, (sptimes[spikecount] - time));
				}
				else {
					sptimes[spikecount] = time;
				}
			}
			++spikecount;			//incremented at the end so it can be used as position in sptimes			
		}
		
		if (spikecount > OFFSET) {	//finding time of threshold crossings for snapshot
			if (fthresh == -1.0 && vout[0] >= STHRESHOLD && v[0] < STHRESHOLD) {
				if (INTERPOLATE) {
					fthresh = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
				}
				else {
					fthresh = time;
				}
				spike.init[0] = -50.0;
				for (i = 1; i < N; ++i) {			//puts initial state variables into spike template
					spike.init[i] = vout[i];
				}
				for (i = 0; i < dsteps; ++i) {
					spike.ibuf[i] = buf[i];			//puts initial buffer into spike template
				}
				spike.bufpos = bufpos;
			}
			else if (fthresh != -1.0 && sndthresh == -1.0 && vout[0] <= STHRESHOLD && v[0] > STHRESHOLD) {
				if (INTERPOLATE) {
					sndthresh = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
				}
				else {
					sndthresh = time;
				}
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
	if (PLONG) {
		makedata(y, xx, nstep, V, "v.data");
		makedata(y, xx, nstep, M, "m.data");
		makedata(y, xx, nstep, H, "h.data");
		makedata(y, xx, nstep, NV, "n.data");
	}
	else {
		printf("\n\nSince PLONG == 0, v-n.data are not being written\n\n");
	}
	dump_(vout);
	
	printf("This simulation counted %d spikes in all.\n", spikecount);
	if (spikecount >= (SAMPLESIZE + OFFSET)) {
		for (i = OFFSET; i < SAMPLESIZE + OFFSET - 1; ++i) {		//calculates differences between spike times to find each period
			sumdiffs += sptimes[i + 1] - sptimes[i];
			spdiffs[i - OFFSET] = sptimes[i + 1] - sptimes[i];
		}
		//~ printdarr(spdiffs, SAMPLESIZE - 1);
		printperiod(spdiffs, SAMPLESIZE - 1, "period.data");
		normalperiod = sumdiffs / (SAMPLESIZE - 1);
		psteps = (int)round(normalperiod / STEPSIZE);
		periodsd = calculateSD(spdiffs, SAMPLESIZE - 1);
		printf("The average unperturbed period is %f, which is approximately %d steps.\n", normalperiod, psteps);
		printf("The standard deviation is %f.\n", periodsd);
	}
	else {
		printf("There are not enough spikes to account for sample size and offset or something else has gone wrong.\n");
		printf("Killing because it hasn't passed the test of spikecount >= (SAMPLESIZE + OFFSET).\n");
		printf("v.data-n.data as well as end.data were still written, but no trace or prc processes occured.\n");
		return 0;
	}
	printf("fthresh = %f and sndthresh = %f\n", fthresh, sndthresh);
	
	snapshot(y, xx, nstep, V, fthresh, sndthresh, &spike);
	printemp(&spike);
	/*
	int prcsteps;
	extern int prcmode;
	double phase;
	prcsteps = psteps * 5;
	double targphase;
	int targstep;
	int pertpos;
	int flag;		//perturbation happened
	double flag1;		//first period complete after perturbation
	double flag2;		//second period complete after perturbation
	//~ int flag3;		//third period complete after perturbation, currently unused
	
	if (DO_TRACE) {
		if (!UNPERT) {
			prcmode = True;
		}
		else {
			prcmode = False;
			prcsteps = psteps;
		}
		Phipair trace;
		targphase = TPHASE;
		targstep = (PRCSKIP) ? (psteps * (1 + targphase)) : (psteps * targphase);
		trace.phase = targphase;
		flag1 = 0.0;
		flag2 = 0.0;
		
		for (i = 0; i < N; ++i) {
			dv[i] = 0.0;
			v[i] = 0.0;
			vout[i] = 0.0;
			current[i] = 0.0;
		}
		
		printf("Attempting to use template\n");
		printemp(&spike);
		copyab(spike.init, v, N);
		//~ printdarr(v, N);
		copyab(spike.ibuf, buf, (int)(DELAY / STEPSIZE));
		time = 0;
		bufpos = spike.bufpos;
		del = &buf[bufpos]; //moves the pointer one step ahead in the buffer
		pert = spike.volts;
		pertpos = 0;
		//~ printf("pert points to %p which holds %f\n", pert, *pert);
		
		derivs(time, v, dv, del);	//does running this one effect the phase/ perturbation? I don't think so but I'm not sure.
	
		rk4(v, dv, N, time, STEPSIZE, vout, del);
		
		xx[0] = time;
		for (i = 0; i < N; i++) {
			v[i] = vout[i]; 
			y[0][i] = v[i];
		}
		
		//~ bufpos = spike.bufpos;
		for (k = 0; k < prcsteps; k++) {
			
			del = &buf[bufpos]; //moves the pointer one step ahead in the buffer
			derivs(time, v, dv, del);
			rk4(v, dv, N, time, STEPSIZE, vout, del);
			*del = vout[0];
				
			time += STEPSIZE;
			xx[k + 1] = time;
		
			if (bufpos < dsteps - 1) {	//increments bufpos within the buffer each step
				bufpos++;
			}
			else {
				bufpos = 0;
			}
			if (!UNPERT) {
				if (k == targstep) {	//activates perturbation mode if on correct step, allows derivs() to start using the "perturbation synapse" (a [pre-recorded stimulus of the same identical neuron)
					//~ printf("%f\n", time);
					pertmode = True;
					flag = True;
				}
			}
			if (pertmode) {
				if (pertpos < spike.steps - 1) {
					pertpos++;
					pert = &spike.volts[pertpos];
				}
				else {
					pertmode = False;
				}
			}
			if (flag && vout[0] >= THRESHOLD && v[0] < THRESHOLD && flag1 == 0.0) {
				
				if (INTERPOLATE) {
					flag1 = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
					trace.fphi1 = (PRCSKIP) ? ((((double)(flag1) - (double)(normalperiod)) - (double)(normalperiod)) / (double)(normalperiod)) : ((double)(flag1) - (double)(normalperiod))/ (double)(normalperiod);
				}
				else {
					flag1 = k;
					trace.fphi1 = (PRCSKIP) ? ((((double)(k) - (double)(psteps)) - (double)(psteps)) / (double)(psteps)) : ((double)(k) - (double)(psteps))/ (double)(psteps);
				}
				printf("The trace was done at phase %f. The step targeted was %d (time = %f), and the total number of steps in the unperturbed period was %d.\n", TPHASE, targstep, time, psteps);
				printf("f(phi)1 is %f.\n", trace.fphi1);
			}
			else if (flag && vout[0] >= THRESHOLD && v[0] < THRESHOLD && flag1 != 0.0 && flag2 == 0.0) {
				if (INTERPOLATE) {
					flag2 = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
					trace.fphi2 = (PRCSKIP) ? (((((double)(flag2) - (double)(flag1)) - (double)(normalperiod)) - (double)(flag1)) / (double)(normalperiod)) : ((double)(flag2) - (double)(flag1) - (double)(normalperiod))/ (double)(normalperiod);
				}
				else {
					flag2 = k;
					trace.fphi2 = ((double)(k) - (double)(flag1))/ (double)(psteps);
					//~ printf("k %d ptime %f targstep %d prc[prcpos].phase %f time %f prc[prcpos].fphi1 %f flag1 %d flag2 %d prc[prcpos].fphi2 %f\n", k, ptime, targstep, prc[prcpos].phase, time, prc[prcpos].fphi1, flag1, flag2, prc[prcpos].fphi2);
				}
			}
			for (i = 0; i < N; i++) {
				v[i] = vout[i];
				y[k + 1][i] = v[i];
			}
		}
		
		makedata(y, xx, prcsteps, V, "tracev.data");	
		makedata(y, xx, prcsteps, M, "tracem.data");
		makedata(y, xx, prcsteps, H, "traceh.data");
		makedata(y, xx, prcsteps, NV, "tracen.data");
	}
	if (DO_PRC) {
		int j;
		Phipair prc[INTERVAL + 1];
		int prcpos = 0;
		
		double ptime;

		for (j = 0; j < INTERVAL + 1; ++j) {
			targphase = ((double)j) * (1.0 / INTERVAL);
			targstep = (PRCSKIP) ? (psteps * (1 + targphase)) : (psteps * targphase);
			
			
			flag = False;
			flag1 = 0.0;
			flag2 = 0.0;
			//~ flag3 = False;
			
			for (i = 0; i < N; ++i) {
				dv[i] = 0.0;
				v[i] = 0.0;
				vout[i] = 0.0;
				current[i] = 0.0;
			}
			
			//~ printf("Attempting to use template\n");
			//~ printemp(&spike);
			copyab(spike.init, v, N);
			//~ printdarr(v, N);
			copyab(spike.ibuf, buf, (int)(DELAY / STEPSIZE));
			time = 0;
			bufpos = spike.bufpos;
			del = &buf[bufpos]; //moves the pointer one step ahead in the buffer
			
			pert = spike.volts;
			pertpos = 0;
			//~ printf("Phase %f Step %d\n", targphase, targstep);
			
			derivs(time, v, dv, del);	//does running this one effect the phase/ perturbation? I don't think so but I'm not sure.
	
			rk4(v, dv, N, time, STEPSIZE, vout, del);
			
			xx[0] = time;
			for (i = 0; i < N; i++) {
				v[i] = vout[i]; 
				y[0][i] = v[i];
			}

			for (k = 0; k < prcsteps; k++) {
				
				del = &buf[bufpos]; //moves the pointer one step ahead in the buffer
				derivs(time, v, dv, del);
				rk4(v, dv, N, time, STEPSIZE, vout, del);
				*del = vout[0];
					
				time += STEPSIZE;
				xx[k + 1] = time;
			
				if (bufpos < dsteps - 1) {	//increments bufpos within the buffer each step
					bufpos++;
				}
				else {
					bufpos = 0;
				}
				
				if (k == targstep) {
					pertmode = True;
					flag = True;
					ptime = time;
				}
				if (pertmode) {
					if (pertpos < spike.steps - 1) {
						pertpos++;
						pert = &spike.volts[pertpos];
					}
					else {
						pertmode = False;
					}
				}
				if (flag && vout[0] >= THRESHOLD && v[0] < THRESHOLD && flag1 == 0.0) {
					
					//~ if (INTERPOLATE) {
						//~ flag1 = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
						//~ trace.fphi1 = (PRCSKIP) ? ((((double)(flag1) - (double)(normalperiod)) - (double)(normalperiod)) / (double)(normalperiod)) : ((double)(flag1) - (double)(normalperiod))/ (double)(normalperiod);
					//~ }
					
					if (INTERPOLATE) {
						flag1 = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
						prc[prcpos].fphi1 = (PRCSKIP) ? ((((double)(flag1) - (double)(normalperiod)) - (double)(normalperiod)) / (double)(normalperiod)) : ((double)(flag1) - (double)(normalperiod)) / (double)(normalperiod);
						prc[prcpos].phase = targphase;
						//~ printf("prc[prcpos].fphi1 %f\n\n", prc[prcpos].fphi1);
					}
					else {
						flag1 = k;
						prc[prcpos].fphi1 = (PRCSKIP) ? ((((double)(k) - (double)(psteps)) - (double)(psteps)) / (double)(psteps)) : ((double)(k) - (double)(psteps))/ (double)(psteps);
						prc[prcpos].phase = targphase;
						//~ ((double)targstep - (double)psteps) / (k - psteps);
					}
					//~ prcpos++;
				}
				else if (flag && vout[0] >= THRESHOLD && v[0] < THRESHOLD && flag1 != 0.0 && flag2 == 0.0) {
					if (INTERPOLATE) {
						flag2 = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
						prc[prcpos].fphi2 = ((double)(flag2) - (double)(flag1) - (double)(normalperiod)) / (double)(normalperiod);
						printf("k %d prcpos %d ptime %f targstep %d prc[prcpos].phase %f time %f prc[prcpos].fphi1 %f flag1 %f flag2 %f prc[prcpos].fphi2 %f\n", k, prcpos, ptime, targstep, prc[prcpos].phase, time, prc[prcpos].fphi1, flag1, flag2, prc[prcpos].fphi2);
						prcpos++;
					}
					else {
						flag2 = k;
						prc[prcpos].fphi2 = ((double)(k) - (double)(flag1))/ (double)(psteps);
						printf("k %d prcpos %d ptime %f targstep %d prc[prcpos].phase %f time %f prc[prcpos].fphi1 %f flag1 %f flag2 %f prc[prcpos].fphi2 %f\n", k, prcpos, ptime, targstep, prc[prcpos].phase, time, prc[prcpos].fphi1, flag1, flag2, prc[prcpos].fphi2);
						prcpos++;
					}
				}
				for (i = 0; i < N; i++) {
					v[i] = vout[i];
					y[k + 1][i] = v[i];
				}
			}
		}
		printphi(prc, INTERVAL, 1, "prc1.data");
		printphi(prc, INTERVAL, 2, "prc2.data");
	}
	*/
	
	if (DO_PRC) {
		prc(spike, INTERVAL, normalperiod);
	}
	if (DO_TRACE) {
		Phipair test;
		test.phase = -1.0;
		printf("test.fphi1 = %f\n", test.fphi1);
		pertsim(normalperiod, spike, &test, 1, "unpert");
		Phipair trace;
		trace.phase = 0.0;
		pertsim(normalperiod, spike, &trace, 1, "0");
		trace.phase = 0.5;
		//~ pertsim(normalperiod, spike, &trace, 1, "0.5");
		//~ trace.phase = 0.6;
		pertsim(normalperiod, spike, &trace, 1, "0.6");
		
		test.fphi1 = 0.0;
		test.fphi2 = 0.0;
		printf("test.fphi1 = %f\n", test.fphi1);	
		printf("trace.fphi1 = %f\n", trace.fphi1);	
	}
	for (i = 0; i < (nstep + 1); i++) {		
		free(y[i]);
	}
	free(xx);	
	
	free(spike.volts);
	
	return 0;
}
