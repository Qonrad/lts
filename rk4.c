//Runge-Kutta Differential Equation Solver

//#include "nrutil.h"
//#include "nlts.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define NN 2
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
#define USE_I_APP 0
#define I_APP_START 500
#define I_APP_END 501
#define USE_LOWPROPOFOL 1	//obviously low and high propofol can't be used together, if both are 1, then lowpropofol is used
#define USE_HIGHPROPOFOL 0
#define PROPOFOL_START 300
#define PROPOFOL_END 100000
#define LOWPROP_GSYN 0.25
#define LOWPROP_TAU 10
#define HIGHPROP_GSYN 0.5
#define HIGHPROP_TAU 20
#define STARTTIME 0
#define ENDTIME 4700
#define STEPSIZE 0.01
#define DELAY 1.0 			//delay must evenly divide stepsize, and it is only used if it is >= stepsize
#define THRESHOLD -50.0		//the voltage at which it counts a spike has occured, used to measure both nonperturbed and perturbed period for PRC
#define STHRESHOLD -50.0	//threshold used to measure just the spike, not the period between spikes
#define SAMPLESIZE 5 		//number of spikes that are averaged together to give unperturbed period
#define OFFSET 20			//number of spikes that are skipped to allow the simulation to "cool down" before it starts measuring the period
#define POPULATION 20		//number of neurons in the whole population
#define MYCLUSTER 10		//number of neurons in the simulated neuron's population
#define DO_PRC 0			//toggle for prc
#define DO_TRACE 0			//toggles doing trace for a single (or multiple phase perturbations) but each is recorded individually
#define TPHASE 0.985
#define INTERVAL 100			//number of intervals prc analysis will be done on
#define True 1
#define False 0
#define INTERPOLATE 1
#define PLONG 1
#define SELF 0

double current[C];	//external current variable, similar to how Canavier did it
static double *del;
static int prcmode;
static int pertmode;
double *pert;
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

typedef struct Phipairs {
	double phase;
	double fphi1;
	double fphi2;
} Phipair;

typedef struct Templates {
	int steps;
	double *volts;
	double init[N * NN];
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

void derivs(double time, double *y, double *dydx, double *oldv, double* weight[NN]) { 
	double iapp, synsum, subsum;
	double *pert;
	int i, j;
	
	
	if (USE_I_APP && !(prcmode)) {
		iapp = (time < I_APP_START || time > I_APP_END) ? I_APP : I_APP_STEP;	
	}
	
	else {
		iapp = I_APP;
	}
	
	
	//~ gsyn *=  weight[i][pre];
	for (i = 0; i < NN; i++) {
		//cell is self-connected and the one that is measured for the PRC.
		// (((y[V + (N * i)] - (-54)) / 4) < 10e-6) ? (0.32 * 4.0) :
		current[I_NA] = G_NA * y[H + (N * i)] * pow(y[M + (N * i)], 3.0) * (y[V + (N * i)] - E_NA);
		current[I_K] =  G_K * pow(y[NV + (N * i)], 4.0) * (y[V + (N * i)] - E_K);
		current[I_M] =  G_M * y[MN + (N * i)] * (y[V + (N * i)] - E_K);
		current[I_L] =  G_L * (y[V + (N * i)] - E_L);
		if (NN == 1) {
			current[I_S] =  gsyn * ((y[S + (N * i)] * (MYCLUSTER - 1))+ (y[P + (N * i)] * (POPULATION - MYCLUSTER))) * (y[V + (N * i)] - E_SYN);
		}
		else {
			//???????????vvvvvvv
			current[I_S] = gsyn * y[S + (N * i)] * (y[V + (N * i)] - E_SYN);
		}
		
		dydx[V + (N * i)] = (iapp - current[I_NA] - current[I_K] - current[I_M] - current[I_L] - current[I_S]) / CM;
		dydx[M + (N * i)] =  (( fabs(((y[V + (N * i)] + 54)) / 4) < 10e-6) ? (0.32 * 4.0) : ( f(y[V + (N * i)], 0.32, -54, 4.0) )) * (1.0 - y[M + (N * i)]) - ((fabs(((y[V + (N * i)] + 27)) / 5) < 10e-6) ? (-0.28 * -5) : ( f(y[V + (N * i)], -0.28, -27, -5.0) )) * y[M + (N * i)];
		//0.32 * (y[V + (N * i)] + 54.0) / (1.0 - exp(-(y[V + (N * i)] + 54.0) / 4.0)) * (1.0 - y[M + (N * i)]) - 0.28 * (y[V + (N * i)] + 27.0) / (exp((y[V + (N * i)] + 27.0) / 5.0) - 1.0) * y[M + (N * i)];   
		dydx[H + (N * i)] = 0.128 * exp(-(y[V + (N * i)] + 50.0) / 18.0) * (1.0 - y[H + (N * i)]) - 4.0 / (1.0 + exp(-(y[V + (N * i)] + 27.0) / 5.0)) * y[H + (N * i)];   
		dydx[NV + (N * i)] = ((fabs(((y[V + (N * i)] + 52)) / 5) < 10e-6) ? (0.032 * 5) : (f(y[V + (N * i)], 0.032, -52, 5.0) )) * (1.0 - y[NV + (N * i)]) - 0.5 * exp(-(y[V + (N * i)] + 57.0) / 40.0) * y[NV + (N * i)];
		//0.032 * (y[V + (N * i)] + 52.0) / (1.0 - exp(-(y[V + (N * i)] + 52.0) / 5.0)) * (1.0 - y[NV + (N * i)]) - 0.5 * exp(-(y[V + (N * i)] + 57.0) / 40.0) * y[NV + (N * i)];  
		dydx[MN + (N * i)] = ((fabs(((y[V + (N * i)] + 30)) / 9) < 10e-6) ? (3.209 * 0.0001 * 9.0) : (f(y[V + (N * i)], (3.209 * 0.0001), -30, 9.0) )) * (1.0 - y[MN + (N * i)]) + (((fabs(((y[V + (N * i)] + 30)) / 9) < 10e-6) ? (3.209 * 0.0001 * -9.0) : (f(y[V + (N * i)], (3.209 * 0.0001), -30, -9.0) )) * y[MN + (N * i)]);
		//above has -q in both approximation and formula, reexamine this!!!!
		//3.209 * 0.0001 * ((y[V + (N * i)] + 30.0) / (1.0 - exp(-(y[V + (N * i)] + 30.0) / 9.0)) * (1.0 - y[MN + (N * i)]) + (y[V + (N * i)] + 30.0) / (1.0 - exp((y[V + (N * i)] + 30.0) / 9.0)) * y[MN + (N * i)]); 
		
		//~ fprintf(stderr, "%lf\n", synsum);
		if (DELAY >= STEPSIZE) {
			//newest edit: *oldv should be the 2 * tanh() function of the previous voltage, moved outside of the driver to prevent waste of computational resources
			synsum = 0.0;
			for (j = 0; j < NN; ++j) {
					synsum += oldv[j] * weight[i][j];
			}

			//?????????????????vvvvvvv
			dydx[S + (N * i)] = (synsum * (1 - y[S + (N * i)])) - (y[S + (N * i)] / tau); //uses *oldv which should be del, the delay pointer in the buffer
			//~ dydx[S + (N * i)] = 2 * (1 + tanh(*oldv / 4.0)) * (1 - y[S + (N * i)]) - y[S + (N * i)] / tau; //uses *oldv which should be del, the delay pointer in the buffer
			//~ printf("I exist at time %f.\n", time);
		}
		else {
			dydx[S + (N * i)] = 2 * (1 + tanh(y[V + (N * i)] / 4.0)) * (1 - y[S + (N * i)]) - y[S + (N * i)] / tau;
		}
		if (NN == 1) {
			if (pertmode) {
				dydx[P + (N * i)] = 2 * (1 + tanh(*pert / 4.0)) * (1 - y[P + (N * i)]) - y[P + (N * i)] / tau;	//should probably be Ps instead of Ss
			}
			else {
				dydx[P + (N * i)] = 2 * (1 + tanh((THRESHOLD) / 4.0)) * (1 - y[P + (N * i)]) - y[P + (N * i)] / tau;
			}
		}
		else {
			dydx[P + (N * i)] = 0;
		}
		//2*(1+tanh(y[V1 + N*j + (N * i)]/4.0))*(1-y[S1 + N*j + (N * i)])-y[S1 + N*j]/TAUSYN;  
	}
	return;
}

void scan_(double *Y) {
	FILE *fopen(),*sp;
	int i, j;
	sp = fopen("state.data","r");
/*	for (i = 0; i < (N); ++i) {
		fscanf(sp, "%lf\n", &Y[i]);
		if (NN > 1) {
			for (j = 0; j < NN; ++j) {
				Y[i + (N * j)] = Y[i];
			}
		}
		fprintf(stderr, "%lf\n", Y[i]);
	}
*/ 
	for (i = 0; i < N*NN; ++i)
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
			fprintf(stderr, "%g\n", Y[pos]);
			++pos;
		}
	}	
	fclose(sp);
}

void rk4(double y[], double dydx[], int n, double x, double h, double yout[], double *oldv, double* weight[NN]) {
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

/*
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
	printdarr(temp->init, (N * NN));
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

void pertsim(double normalperiod, Template spike, Phipair *trace, int tracedata, const char* tracename) {
	int i, k, nstep, targstep, flag, pertpos, psteps;
	double time, flag1, flag2;
	double v[N * NN], vout[N * NN], dv[N * NN];			//v = variables (state and current), vout = output variables, dv = derivatives (fstate)
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
		y[i] = (double*) malloc(sizeof(double) * (N * NN));
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
	for (i = 0; i < (N * NN); ++i) {
			dv[i] = 0.0;
			v[i] = 0.0;
			vout[i] = 0.0;
			current[i] = 0.0;
	}
	
	//~ printf("\n\n\n\nAttempting to use template within pertsim!\n");
	//~ printemp(&spike);
	copyab(spike.init, v, (N * NN));
	copyab(spike.ibuf, buf, (int)(DELAY / STEPSIZE));
	bufpos = spike.bufpos;
	del = &buf[bufpos]; //moves the pointer to the correct initial position in the buffer, unnecessary i believe
	pert = spike.volts;
	pertpos = 0;
	//~ printf("Template successfully initiated!\n");
	
	//~ derivs(time, v, dv, del);	//does running this one effect the phase/ perturbation? I don't think so but I'm not sure.
	//~ rk4(v, dv, (N * NN), time, STEPSIZE, vout, del);
	
	xx[0] = time;
	for (i = 0; i < (N * NN); i++) {
		//~ v[i] = vout[i]; 
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
		rk4(v, dv, (N * NN), time, STEPSIZE, vout, del);
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
		for (i = 0; i < (N * NN); i++) {
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

void makeunpert(double** y, double *xx, int normalperiod, int startstep, int realtime, int full, const char *filename) {
	int i, psteps, x;
	psteps = (int)round(normalperiod / STEPSIZE);
	FILE *fopen(),*fp;
	fp = fopen(filename, "w");
	if (realtime) {
		if (!full) {
			for (i = startstep; i < (startstep + (5 * psteps)); i++) {
				fprintf(fp, "%f %f\n", xx[i], y[i][V]);
			}
		}
		else {
			for (i = startstep; i < (startstep + (5 * psteps)); i++) {
				fprintf(fp, "%f %f %f %f %f %f %f %f\n", xx[i], y[i][0], y[i][1], y[i][2], y[i][3], y[i][4], y[i][5], y[i][6]);
			}
		}
	}
	else {
		if (full) {
			for (i = startstep; i < (startstep + (5 * psteps)); i++) {
				fprintf(fp, "%f %f %f %f %f %f %f %f\n", xx[x], y[i][0], y[i][1], y[i][2], y[i][3], y[i][4], y[i][5], y[i][6]);
				x++;
			}
		}
		else {
			for (i = startstep; i < (startstep + (5 * psteps)); i++) {
				fprintf(fp, "%f %f\n", xx[x], y[i][V]);
				x++;
			}
		}
	}
	fclose(fp);
}
*/
int main() {
	gsyn = (G_SYN);
	tau = TAUSYN;
	//Variables to do with simulation
	long int nstep = (ENDTIME - STARTTIME) / STEPSIZE;	// This assumes the endtime-starttime is evenly divisible by the stepsize, which should always be true I think
	printf("The initial simulation will contain %d steps.\n", nstep);
	int i, j, k, pre;
	double time;
	double v[N * NN], vout[N * NN], dv[N * NN];					//v = variables (state and current), vout = output variables, dv = derivatives (fstate)
	double **y, *xx; 						//results variables, y[1..N][1..NSTEP+1], xx[1..NSTEP+1]
	extern double current[];				//external variable declaration
	double* weight[NN];
	for (i = 0; i < NN; i++) {
		weight[i] = (double*) malloc(sizeof(double) * NN);
		for (pre = 0; pre < NN; ++pre) {
			weight[i][pre] = 1.0 / ((double)(NN));
		}
		if (!SELF) {
			weight[i][i] = 0.;
		}
	}
	fprintf(stderr, "test1\n");
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
	
	//Variables for storing spike snapshot
	double fthresh = -1.0;					//time at which voltage crosses threshold on the way up
	double sndthresh = -1.0;				//time at which voltage crosses threshold on the way down
	Template spike;

	int startstep;
	fprintf(stderr, "test2\n");
	
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
	
	//Variables to do with conducting prc
	extern int pertmode;
	extern double *pert;
	
	
	
	time = STARTTIME;
	scan_(v);				//scanning in initial variables (state variables only) 
	fprintf(stderr, "test3\n");
	printdarr(v, (N * NN));
	for (i = 0; i < (dsteps); ++i) {//sets every double in buffer(s) to be equal to the steady state (initial) voltage that was just scanned in
		for (j = 0; j < NN; ++j) {
			if (DELAY >= STEPSIZE) { 
				buf[i][j] = 2 * (1 + tanh(v[V + (N * j)] / 4.0));
			}
			else {
				buf[i][j] = v[V + (N * j)];
			}
		}
	}
	fprintf(stderr, "test4\n");
	derivs(time, v, dv, del, weight);
	printdarr(dv, (N * NN));
	rk4(v, dv, (N * NN), time, STEPSIZE, vout, del, weight);
	
	
	xx[0] = STARTTIME;		
	for (i = 0; i < (N * NN); i++) {
		v[i] = vout[i]; 
		y[0][i] = v[i];
		
	}
	fprintf(stderr, "test5\n");
	printdarr(v, (N * NN));
	for (k = 0; k < nstep; k++) {
		
		if (USE_LOWPROPOFOL) {	//changes gsyn to match the correct level of propofol in the simulation
			gsyn = (time > PROPOFOL_START || time < PROPOFOL_END || prcmode) ? LOWPROP_GSYN : G_SYN;
			tau = (time > PROPOFOL_START || time < PROPOFOL_END || prcmode) ? LOWPROP_TAU : TAUSYN; 
		}
		else if (USE_HIGHPROPOFOL) {
			gsyn = (time > PROPOFOL_START || time < PROPOFOL_END || prcmode) ? HIGHPROP_GSYN : G_SYN;
			tau = (time > PROPOFOL_START || time < PROPOFOL_END || prcmode) ? HIGHPROP_TAU : TAUSYN;
		}
		else {	
			gsyn = (G_SYN);
			tau = TAUSYN;
		}
		
		del = buf[bufpos]; //moves the pointer one step ahead in the buffer
		//~ printdarr(v, (N * NN));
		//~ fprintf(stderr, "test5a k = %d\n", k);
		derivs(time, v, dv, del, weight);										//actually does the step
		//~ printdarr(v, (N * NN));
		//~ fprintf(stderr, "test5b k = %d\n", k);
		rk4(v, dv, (N * NN), time, STEPSIZE, vout, del, weight);		//actually dose the step
		//~ printdarr(v, (N * NN));
		//~ fprintf(stderr, "test5c k = %d\n", k);
		for (j = 0; j < NN; ++j) {								//IMPORTANT: this calculates the f(vj) outside of derivs and places that into the buffer, preventing a large waste of computational time
			del[j] = 2 * (1 + tanh(vout[V + (N * j)] / 4.0));		//dereferences the delay pointer, and puts the previous f(vj) in the same spot
		}
		/*
		if (DO_TRACE || DO_PRC) {
			if (vout[0] >= THRESHOLD && v[0] < THRESHOLD) {
				if (spikecount < (SAMPLESIZE + OFFSET)) {
					if (INTERPOLATE) {
						sptimes[spikecount] = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
						vout[0] = -50.0; //fudging the interpolated value to the most recent v[0] value so that it exactly matches the spike template
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
					printf("Using variables at time %f (interpolated time %f) for start of spike template.\n", time, fthresh);
					printf("Current state variables, except voltage is interpolated to 0.\n");
					printdarr(vout, (N * NN));
					startstep = k + 1;
					
					for (i = 1; i < (N * NN); ++i) {			//puts initial state variables into spike template
						spike.init[i] = vout[i];
					}
					printf("Printing the spike.init array just in case.\n");
					printdarr(spike.init, (N * NN));
					
					printf("About to copy buffer and buffer position to template.\n");
					printf("Buffer array.\n");
					printdarr(buf, dsteps);
					printf("bufpos = %f\n", bufpos);
					
					for (i = 0; i < dsteps; ++i) {
						spike.ibuf[i] = buf[i];			//puts initial buffer into spike template
					}
					spike.bufpos = bufpos;
					
					printf("Done copying. spike.ibuf\n");
					printdarr(spike.ibuf, dsteps);
					printf("spike.bufpos = %f\n", spike.bufpos);
					
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
		
		}		
		*/
		
		time += STEPSIZE;
		xx[k + 1] = time;
		
/*		if (bufpos < dsteps - 1) {	//increments bufpos within the buffer each step
			bufpos++;
		}
		else {
			bufpos = 0;
		}
*/
		bufpos = (++bufpos)%dsteps;
		
		for (i = 0; i < (N * NN); i++) {
			v[i] = vout[i];
			y[k + 1][i] = v[i];
		}
	}
	fprintf(stderr, "test6\n");
	if (PLONG) {
		makedata(y, xx, nstep, V, "v.data");
		makedata(y, xx, nstep, (V + N), "v2.data");
		makedata(y, xx, nstep, M, "m.data");
		makedata(y, xx, nstep, H, "h.data");
		makedata(y, xx, nstep, NV, "n.data");
	}
	else {
		printf("\n\nSince PLONG == 0, v-n.data are not being written\n\n");
	}
	fprintf(stderr, "test7\n");
	printdarr(vout, (N * NN));
	dump_(vout);
/*	
	if (DO_TRACE || DO_PRC) {
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
		makeunpert(y, xx, normalperiod, startstep, 0, 1, "unpertvx.data");
	}
	

	snapshot(y, xx, nstep, V, fthresh, sndthresh, &spike);
	printemp(&spike);

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
		pertsim(normalperiod, spike, &trace, 1, "0.5");
		trace.phase = 0.6;
		pertsim(normalperiod, spike, &trace, 1, "0.6");
		
		test.fphi1 = 0.0;
		test.fphi2 = 0.0;
		printf("test.fphi1 = %f\n", test.fphi1);	
		printf("trace.fphi1 = %f\n", trace.fphi1);	
	}
	* 
*/

	for (i = 0; i < (nstep + 1); i++) {		
		free(y[i]);
	}
	free(xx);	
	
	//~ free(spike.volts);
	
	return 0;
}
