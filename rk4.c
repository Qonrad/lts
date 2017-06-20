//Runge-Kutta Differential Equation Solver

//#include "nrutil.h"
//#include "nlts.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define N 6
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
#define CM 1.00 /*uF/cm2*/
#define I_APP 1.8 /*uA/cm2*/ //Conrad, should be 2.0
#define E_NA  50.0
#define E_K  -100.0
#define E_L  -67.0
#define E_SYN  -80.0
#define G_NA 100.0 /* mS/cm2*/
#define G_K   80.0
#define G_M   2	// Was previously almost always 2, McCarthy seems to have it at 4, gmi
#define G_L   0.1
#define G_SYN  0.5	//McCarthy gi_i baseline = 0.165, low-dose Propofol = 0.25, high-dose Propofol = 0.5
#define TAUSYN 20		//McCarthy taui baseline = 5.0, low-dose Propofol = 10, high-dose Propofol = 20
#define USE_I_APP 1
#define STARTTIME 1
#define ENDTIME 1200	
#define STEPSIZE 0.05

double current[C];	//external current variable, similar to how Canavier did it
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
	printf("[");
	for (i = 0; i < numelems - 1; ++i) {
		printf("%f (%p), ", a[i], &a[i]);
	}
	printf("%f (%p)]\n", a[numelems - 1], &a[numelems - 1]);
}

inline double f(double v, double a, double th, double q) {
	return a * ((v - th) / (1 - exp(-(v - th) / q)));
}

derivs(double time, double *y, double *dydx) { 
	double iapp;
	extern double current[];
	
	if (USE_I_APP) {
		iapp = (time < 300 || time > 330) ? I_APP : 3.2;
	}
	else {
		iapp = I_APP;
	}
	// (((y[V] - (-54)) / 4) < 10e-6) ? (0.32 * 4.0) :
	current[I_NA] = G_NA * y[H] * pow(y[M], 3.0) * (y[V] - E_NA);
	current[I_K] =  G_K * pow(y[NV], 4.0) * (y[V] - E_K);
	current[I_M] =  G_M * y[MN] * (y[V] - E_K);
	current[I_L] =  G_L * (y[V] - E_L);
	current[I_S] =  G_SYN * y[S] * (y[V] - E_SYN);
	
	dydx[V] = (iapp - current[I_NA] - current[I_K] - current[I_M] - current[I_L] - current[I_S]) / CM;
	dydx[M] =  (( fabs(((y[V] + 54)) / 4) < 10e-6) ? (0.32 * 4.0) : ( f(y[V], 0.32, -54, 4.0) )) * (1.0 - y[M]) - ((fabs(((y[V] + 27)) / 5) < 10e-6) ? (-0.28 * -5) : ( f(y[V], -0.28, -27, -5.0) )) * y[M];
	//0.32 * (y[V] + 54.0) / (1.0 - exp(-(y[V] + 54.0) / 4.0)) * (1.0 - y[M]) - 0.28 * (y[V] + 27.0) / (exp((y[V] + 27.0) / 5.0) - 1.0) * y[M];   
	dydx[H] = 0.128 * exp(-(y[V] + 50.0) / 18.0) * (1.0 - y[H]) - 4.0 / (1.0 + exp(-(y[V] + 27.0) / 5.0)) * y[H];   
	dydx[NV] = ((fabs(((y[V] + 52)) / 5) < 10e-6) ? (0.032 * 5) : (f(y[V], 0.032, -52, 5.0) )) * (1.0 - y[NV]) - 0.5 * exp(-(y[V] + 57.0) / 40.0) * y[NV];
	//0.032 * (y[V] + 52.0) / (1.0 - exp(-(y[V] + 52.0) / 5.0)) * (1.0 - y[NV]) - 0.5 * exp(-(y[V] + 57.0) / 40.0) * y[NV];  
	dydx[MN] = ((fabs(((y[V] + 30)) / 9) < 10e-6) ? (3.209 * 0.0001 * 9.0) : (f(y[V], (3.209 * 0.0001), -30, 9.0) )) * (1.0 - y[MN]) + (((fabs(((y[V] + 30)) / 9) < 10e-6) ? (3.209 * 0.0001 * -9.0) : (f(y[V], (3.209 * 0.0001), -30, -9.0) )) * y[MN]);
	//above has -q in both approximation and formula, reexamine this!!!!
	//3.209 * 0.0001 * ((y[V] + 30.0) / (1.0 - exp(-(y[V] + 30.0) / 9.0)) * (1.0 - y[MN]) + (y[V] + 30.0) / (1.0 - exp((y[V] + 30.0) / 9.0)) * y[MN]); 
	
	dydx[S] = 2 * (1 + tanh(y[V] / 4.0)) * (1 - y[S]) - y[S] / TAUSYN; //2*(1+tanh(y[V1 + N*j]/4.0))*(1-y[S1 + N*j])-y[S1 + N*j]/TAUSYN;  
	
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

void rk4(double y[], double dydx[], int n, double x, double h, double yout[]) {
	int i;
	double xh, hh, h6, *dym, *dyt, *yt;
	
	dym = (double*) malloc(n * sizeof(double));
	dyt = (double*) malloc(n * sizeof(double));
	yt = (double*) malloc(n * sizeof(double));
	
	
	hh = h * 0.5;
	h6 = h / 6.0;
	xh = x + hh;
	//~ printf("before 1st step, uses y, hh, dydx, adds to yt time = %f\n", x); //these are all print statements for debugging
	//~ printf("y\n");
	//~ printdarr(y, N);
	//~ printf("yt\n");
	//~ printdarr(yt, N);
	//~ printf("dydx\n");
	//~ printdarr(dydx, N);
	//~ printf("dyt\n");
	//~ printdarr(dyt, N);
	for (i = 0; i < n; i++) {			//first step
		yt[i] = y[i] + hh * dydx[i];
	}
	//~ printf("before 1st step, derivs(xh, yt, dyt); time = %f\n", x);
	//~ printf("y\n");
	//~ printdarr(y, N);
	//~ printf("yt\n");
	//~ printdarr(yt, N);
	//~ printf("dydx\n");
	//~ printdarr(dydx, N);
	//~ printf("dyt\n");
	//~ printdarr(dyt, N);
	derivs(xh, yt, dyt);				//second step
	
	for (i = 0; i < n; i++) {
		yt[i] = y[i] + hh * dyt[i];
	}
	//~ printf("before 3rd step, uses xh, yt, dym time = %f\n", x);
	//~ printf("y\n");
	//~ printdarr(y, N);
	//~ printf("yt\n");
	//~ printdarr(yt, N);
	//~ printf("dydx\n");
	//~ printdarr(dydx, N);
	//~ printf("dyt\n");
	//~ printdarr(dyt, N);
	//~ printf("dym\n");
	//~ printdarr(dym, N);
	derivs(xh, yt, dym);				//third step
	
	for (i = 0; i < n; i++) {
		yt[i] = y[i] + h * dym[i];
		dym[i] += dyt[i];
	}
	//~ printf("before 4th step, derivs(x + h, yt, dyt); time = %f\n", x);
	//~ printf("y\n");
	//~ printdarr(y, N);
	//~ printf("yt\n");
	//~ printdarr(yt, N);
	//~ printf("dydx\n");
	//~ printdarr(dydx, N);
	//~ printf("dyt\n");
	//~ printdarr(dyt, N);
	derivs(x + h, yt, dyt);			//fourth step
	
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
		//~ printf("%f %f\n", xx[i], y[i][var]);
		fprintf(fp, "%f %f\n", xx[i], y[i][var]);
	}
	fclose(fp);
}

int main() {
	int i, k;
	double time;
	double *v, *vout, *dv;	//v = variables (state and current), vout = output variables, dv = derivatives (fstate)
	double **y, *xx; 		//results variables, y[1..N][1..NSTEP+1], xx[1..NSTEP+1]
	int nstep;				//number of steps necessary to reach ENDTIME from STARTTIME at the defined STEPSIZE
	extern double current[];	//external variable declaration
	
	nstep = (ENDTIME - STARTTIME) / STEPSIZE;	// This assumes the entime is evenly divisible by the stepsize, which should always be true I think
/*	y = (double**) malloc(sizeof(double) * N * (nstep + 1) + sizeof(double*)*(nstep+1));
	for(i=0;i<(nstep+1);i++)
		y[i]=(double*)(y
*/

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
	derivs(time, v, dv);
	
	rk4(v, dv, N, time, STEPSIZE, vout);
	//~ printdarr(v, N);
	//~ return 0;
	
	
	xx[0] = STARTTIME;		
	for (i = 0; i < N; i++) {
		v[i] = vout[i]; 
		y[0][i] = v[i];
		
	}
	for (k = 0; k < nstep; k++) {
		derivs(time, v, dv);
		rk4(v, dv, N, time, STEPSIZE, vout);
		time += STEPSIZE;
		//~ printf("time = %f\n", time);
		xx[k + 1] = time;
		printf("%f %f\n", time, vout[0]);
		for (i = 0; i < N; i++) {
			//~ printf("v\n");
			//~ printdarr(v, N);
			//~ printf("vout\n");
			//~ printdarr(vout, N);
			v[i] = vout[i];
			y[k + 1][i] = v[i];
		}
		
	}
	//~ printdarr(xx, nstep);
	//~ for (i = 0; i < (nstep + 1); i++) {		//commented out b/c it's causing errors and I don't know why :(
		//~ free(y[i]);
	//~ }
	
	makedata(y, xx, nstep, V, "v.data");
	makedata(y, xx, nstep, M, "m.data");
	makedata(y, xx, nstep, H, "h.data");
	makedata(y, xx, nstep, NV, "n.data");
	
	
	//~ free(current);
	dump_(vout);
	free(v);
	free(dv);
	free(vout);
	//~ free(y);
	free(xx);
	
	return 0;
}
