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
#define DELAY 0.0 			//delay must evenly divide stepsize, and it is only used if it is >= stepsize
#define THRESHOLD -50.0		//the voltage at which it counts a spike has occured, used to measure both nonperturbed and perturbed period for PRC
#define STHRESHOLD -50.0	//threshold used to measure just the spike, not the period between spikes
#define SAMPLESIZE 5 		//number of spikes that are averaged together to give unperturbed period
#define OFFSET 20			//number of spikes that are skipped to allow the simulation to "cool down" before it starts measuring the period
#define POPULATION 20		//number of neurons in the whole population, should be 20 for accurate representation of mccarthy
#define MYCLUSTER 10			//number of neurons in the simulated neuron's population, should be 10 for accurate representation of mccarthy
#define DO_PRC 0			//toggle for prc
#define DO_TRACE 0			//toggles doing trace for a single (or multiple phase perturbations) but each is recorded individually
#define TPHASE 0.985
#define INTERVAL 100			//number of intervals prc analysis will be done on
#define True 1
#define False 0
#define INTERPOLATE 1
#define PLONG 1
#define G(X,Y) ( (fabs((X)/(Y))<1e-6) ? ((Y)*((X)/(Y)/2. - 1.)) : ((X)/(1. - exp( (X)/ (Y) ))) )
#define F(X,Y) ( (fabs((X)/(Y))<1e-6) ? ((Y)*(1.-(X)/(Y)/2.)) : ((X)/(exp( (X)/ (Y) ) -1)) )

double current[C];	//external current variable, similar to how Canavier did it
static double *del;
static int prcmode;
static int pertmode;
double *pert;

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
	fprintf(stderr, "[");
	for (i = 0; i < numelems - 1; ++i) {
		fprintf(stderr, "%f (%p), ", a[i], &a[i]);
	}
	fprintf(stderr, "%f (%p)]\n", a[numelems - 1], &a[numelems - 1]);
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
	
	current[I_NA] = G_NA * y[H] * pow(y[M], 3.0) * (y[V] - E_NA);
	current[I_K] =  G_K * pow(y[NV], 4.0) * (y[V] - E_K);
	current[I_M] =  G_M * y[MN] * (y[V] - E_K);
	current[I_L] =  G_L * (y[V] - E_L);
	current[I_S] =  gsyn * ((y[S] * (MYCLUSTER - 1))+ (y[P] * (POPULATION - MYCLUSTER))) * (y[V] - E_SYN);
	
	dydx[V] = (iapp - current[I_NA] - current[I_K] - current[I_M] - current[I_L] - current[I_S]) / CM;
	dydx[M] = ((0.32) * G((y[V] + 54), -4.) * (1. - y[M])) - ((0.28) * F((y[V] + 27), 5.) * y[M]);
	dydx[H] = 0.128 * exp(-(y[V] + 50.0) / 18.0) * (1.0 - y[H]) - 4.0 / (1.0 + exp(-(y[V] + 27.0) / 5.0)) * y[H];   
	dydx[NV] = ((0.032) * G((y[V] + 52), -5.) * (1. - y[NV])) - (0.5 * exp(-(y[V] + 57.) / 40.) * y[NV]); 
	dydx[MN] = ((3.209 * 0.0001) * G((y[V] + 30), -9.)  * (1.0 - y[MN])) + ((3.209 * 0.0001) * G((y[V] + 30), 9.) * y[MN]);
	
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
	fprintf(stderr, "This template contains %d steps.\n", temp->steps);
	fprintf(stderr, "Array of Voltages\n");
	printdarr(temp->volts, temp->steps);
	fprintf(stderr, "Initial array of state variables\n");
	printdarr(temp->init, N);
	fprintf(stderr, "Initial array of buffer\n");
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

	int startstep;
	
	
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
					vout[0] = -50.0; //fudging the interpolated value to the most recent v[0] value so that it exactly matches the spike template
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
				printdarr(vout, N);
				startstep = k + 1;
				
				for (i = 1; i < N; ++i) {			//puts initial state variables into spike template
					spike.init[i] = vout[i];
				}
				printf("Printing the spike.init array just in case.\n");
				printdarr(spike.init, N);
				
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
		makefull(y, xx, nstep, "nodelfull.data");
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
		printperiod(spdiffs, SAMPLESIZE - 1, "period.data");
		normalperiod = sumdiffs / (SAMPLESIZE - 1);
		psteps = (int)round(normalperiod / STEPSIZE);
		periodsd = calculateSD(spdiffs, SAMPLESIZE - 1);
		printf("The average unperturbed period is %f, which is approximately %d steps.\n", normalperiod, psteps);
		printf("The standard deviation is %f.\n", periodsd);
	}
	else {
		fprintf(stderr, "There are not enough spikes to account for sample size and offset or something else has gone wrong.\n");
		fprintf(stderr, "Killing because it hasn't passed the test of spikecount >= (SAMPLESIZE + OFFSET).\n");
		fprintf(stderr, "v.data-n.data as well as end.data were still written, but no trace or prc processes occured.\n");
		return 0;
	}
	fprintf(stderr, "fthresh = %f and sndthresh = %f\n", fthresh, sndthresh);
	
	
	//~ makeunpert(y, xx, normalperiod, startstep, 0, 1, "unpertvx.data");
	
	
	
	snapshot(y, xx, nstep, V, fthresh, sndthresh, &spike);
	printemp(&spike);
	
	if (DO_PRC) {
		prc(spike, INTERVAL, normalperiod);
	}
	if (DO_TRACE) {
		Phipair test;
		test.phase = -1.0;
		fprintf(stderr, "test.fphi1 = %f\n", test.fphi1);
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
		fprintf(stderr, "test.fphi1 = %f\n", test.fphi1);	
		fprintf(stderr, "trace.fphi1 = %f\n", trace.fphi1);	
	}
	for (i = 0; i < (nstep + 1); i++) {		
		free(y[i]);
	}
	free(xx);	
	
	free(spike.volts);
	
	return 0;
}
