<<<<<<< HEAD
=======
//Runge-Kutta Differential Equation Solver

//#include "nrutil.h"
//#include "nlts.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
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
<<<<<<< HEAD
#define I_APP 2.0 /*uA/cm2*/ //was 1.81 for the prc stuff, experimenting with changing this, risky
=======
#define I_APP 1.81 /*uA/cm2*/ 
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
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
<<<<<<< HEAD
#define TAUSYN 5			//McCarthy taui baseline = 5.0, low-dose Propofol = 10, high-dose Propofol = 20
#define USE_I_APP 0			//really should be called "USE_IAPP_STEP"
#define I_APP_HET 0
#define I_APP_NEURONS 0		//if I_APP enabled, directs the I_APP to affect neurons 0-I_APP_NEURONS 
#define I_APP_START 500
#define I_APP_END 501
#define LOWPROP_GSYN 0.25 //should be 0.25, divided it by 20 instead of using DIVNN to exactly match Carmen's code
=======
#define TAUSYN 10			//McCarthy taui baseline = 5.0, low-dose Propofol = 10, high-dose Propofol = 20
#define USE_I_APP 1
#define I_APP_START 500
#define I_APP_END 501
#define USE_LOWPROPOFOL 1	//obviously low and high propofol can't be used together, if both are 1, then lowpropofol is used
#define USE_HIGHPROPOFOL 0
#define PROPOFOL_START 300
#define PROPOFOL_END 50000
#define LOWPROP_GSYN 0.25
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
#define LOWPROP_TAU 10
#define HIGHPROP_GSYN 0.5
#define HIGHPROP_TAU 20
#define STARTTIME 0
<<<<<<< HEAD
#define DO_TRACE 0			//toggles doing trace for a single (or multiple phase perturbations) but each is recorded individually
#define THRESHOLD -50.0		//the voltage at which it counts a spike has occured, used to measure both nonperturbed and perturbed period for PRC
#define STHRESHOLD -50.0	//threshold used to measure just the spike, not the period between spikes
#define SAMPLESIZE 50 		//number of spikes that are averaged together to give unperturbed period
#define OFFSET 10			//number of spikes that are skipped to allow the simulation to "cool down" before it starts measuring the period
#define True 1
#define False 0
#define INTERPOLATE 1
#define FULLNAME "allvolts.data"
#define DBIT 0
#define G(X,Y) ( (fabs((X)/(Y))<1e-6) ? ((Y)*((X)/(Y)/2. - 1.)) : ((X)/(1. - exp( (X)/ (Y) ))) )
#define F(X,Y) ( (fabs((X)/(Y))<1e-6) ? ((Y)*(1.-(X)/(Y)/2.)) : ((X)/(exp( (X)/ (Y) ) -1)) )
//#define PERTENDTIME 5000	//separate endtime for prc stuff in order to differentiate it from main simulation (seems like a relic of the past?)

//Runge-Kutta Differential Equation Solver, abc

//Attempting to make a unified solution to combine 2neur and fixtemplate branches
//Goal = Be able to put in a single set of parameters and run the code once to produce both a PRC and the true simulations
//Switching between branches as I have it set up currently is tedious and prone to errors. Hopefully this won't be too difficult and will make things much easier.

//Command for compiling and running on macbook
//gcc-7 rk4.c && ./a.out && python2 all-view2.py 20 full.data
//latest update
//gcc-7 -largp rk4.c -o lts && ./lts state.data --prc 100

//Command for compiling and running on beowolf
//gcc rk4.c -lm && ./a.out && python2 all-view2.py 20 full.data
=======
#define ENDTIME 50000
#define STEPSIZE 0.01
#define DELAY 0.0 			//delay must evenly divide stepsize, and it is only used if it is >= stepsize
#define THRESHOLD -50.0		//the voltage at which it counts a spike has occured, used to measure both nonperturbed and perturbed period for PRC
#define STHRESHOLD -50.0	//threshold used to measure just the spike, not the period between spikes
#define SAMPLESIZE 400 		//number of spikes that are averaged together to give unperturbed period
#define OFFSET 200			//number of spikes that are skipped to allow the simulation to "cool down" before it starts measuring the period
#define POPULATION 20		//number of neurons in the whole population
#define MYCLUSTER 10		//number of neurons in the simulated neuron's population
#define DO_PRC 1			//toggle for prc
#define DO_TRACE 0			//toggles doing trace for a single 
#define TPHASE 0.995
#define INTERVAL 200			//number of intervals prc analysis will be done on
#define True 1
#define False 0
#define PRCSKIP 0
#define INTERPOLATE 1
#define UNPERT 1
#define PLONG 0
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0

/*Changes to make
 * randomize inputfile
 * properly able to cat args.txt into parameters to perfectly recreate exact running of specific commit
 * diff viewer
 note: maybe its the divnn that's causing the difference from lowhigh to now, also population & mycluster?
 number of spikes difference with old.c is still sort of mysterious
 * 
 * 
 fixed major bug, now main task is to get those state variables carmen wanted out of the code
 */



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <argp.h>

//prcfuncs
double current[C];	//external current variable, similar to how Canavier did it
static double *del;
<<<<<<< HEAD
double gsyn, tau;
int nn, population;
double *iapps;
static int prcmode;
static int pertmode;
double *pert;
double pweight = 1;
int self_connection = 0;


const char *argp_program_version =
  "lts simulation 1.0";
const char *argp_program_bug_address =
  "<cleonik@tulane.edu>";

/* Program documentation. */
static char doc[] =
  "an lts simulation programmed by Conrad Leonik\n";

/* A description of the arguments we accept. */
static char args_doc[] = "INPUTFILE NUMBER_OF_NUERONS";

/* The options we understand. */
static struct argp_option options[] = {
	{"prc",		'p', "MYCLUSTER",	0, "Run PRC. Argument is size of cluster." },
	{"etime",	'e', "ENDTIME",		0, "ENDTIME of main simulation. Default is 5000."},
	{"lowprop", 'l', "LOW_RANGE",	0, "Use Low-dose Propofol. Mandatory argument is time range eg: 0-5000"},
	{"highprop",'h', "HIGH_RANGE",	0, "Use High-dose Propofol. Mandatory arguement is time range eg: 0-5000"},
	{"delay",	'd', "DELAY",		0, "synaptic delay of lts neurons. default is 0, must evenly divide stepsize"},
	{"stepsize",'s', "STEPSIZE",	0, "size in ms of each step of the simulation, must be evenly divided by the delay if there is one"},
	{"commit",	'c', "toggle",		OPTION_ARG_OPTIONAL, "option to commit the data and code changes when done running"},
	{"verbose", 'v', "toggle",		OPTION_ARG_OPTIONAL, "toggles printing out extra data, only prints spike voltages by default (and prc if that's enabled"},
	{"graph",	'g', "toggle",		OPTION_ARG_OPTIONAL, "toggles the graph option"},
	{"divnn",	'z', "toggle",		OPTION_ARG_OPTIONAL, "toggles dividing by NN, it DOES NOT by default even though McCarthy does"},
	{"interval", 'i', "INTERVAL",	0, "can set interval manually, default is 100"},	 
	{ 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
  char *args[2];                /* arg1 & arg2 */
  int prc;
  int interval;
  int clustersize;
  int etime;
  int lowprop;
  int lowstart;
  int lowend;
  int highprop;
  int highstart;
  int highend;
  int commit;
  int verbose;
  double delay;
  double stepsize;
  int graph;
  int divnn;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
	struct arguments *arguments = state->input;

	switch (key) {
		case 'd':
			arguments->delay = atof(arg);
			break;
		case 'p':
			arguments->prc = 1;
			arguments->clustersize = atoi(arg);
			break;
		case 'e':
			arguments->etime = atoi(arg);
			break;
		case 'h':
			arguments->highprop = 1;
			range_parser(&(arguments->highstart), &(arguments->highend), arg);
			break;
		case 'c':
			arguments->commit = 1;
			break;
		case 's':
			arguments->stepsize = atof(arg);
			break;
		case 'v':
			arguments->verbose = 1;
			break;
		case 'g':
			arguments->graph = 1;
			break;
		case 'l':
			arguments->lowprop = 1;
			range_parser(&(arguments->lowstart), &(arguments->lowend), arg);
			break;
		case 'z':
			arguments->divnn = 1;
			break;
		case 'i':
			arguments->interval = atof(arg);
			break;
		case ARGP_KEY_ARG:
			if (state->arg_num >= 2)
			/* Too many arguments. */
			argp_usage (state);

			arguments->args[state->arg_num] = arg;
			break;

		case ARGP_KEY_END:
			if (state->arg_num < 1)
				/* Not enough arguments. */
				argp_usage (state);
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
  return 0;
=======
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
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
}

<<<<<<< HEAD
/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };
struct arguments arguments;

=======
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
typedef struct Phipairs {
	double phase;
	double fphi1;
	double fphi2;
} Phipair;

typedef struct Templates {
	int steps;
	double *volts;
	double init[N];
<<<<<<< HEAD
	int bufpos;
	double *ibuf;
} Template;

=======
	double ibuf[(int)(DELAY / STEPSIZE)];
	int bufpos;
} Template;

typedef struct TCross {
	int poststep;
	double itime;
} Cross;

//~ void interpolate(Cross* c, prev, postv) {
	//~ 
//~ }

// Conrad code to print an array for debugging
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
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
	double iapp, gsyn, tau;
<<<<<<< HEAD
	int i, j;
	extern double *pert;
	
	
	if (arguments.highprop) {	//this should make it such that highprop will override lowprop
		gsyn = ((time > arguments.highstart && time < arguments.highend)) ? HIGHPROP_GSYN : G_SYN;
		tau = ((time > arguments.highstart && time < arguments.highend)) ? HIGHPROP_TAU : TAUSYN;
	}
	else if (arguments.lowprop) { //is the else truly necessary here?
		gsyn = ((time > arguments.lowstart && time < arguments.lowend)) ? LOWPROP_GSYN : G_SYN;
		tau = ((time > arguments.lowstart && time < arguments.lowend)) ? LOWPROP_TAU : TAUSYN; 
	}
	else {	
		gsyn = G_SYN;
=======
	extern double *pert;
	
	//~ if (pertmode) {
		//~ printf("%f %f\n", time, *pert);
		//~ printf("pert points to %p", pert);
	//~ }
	
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
	}
	else {	
		gsyn = (G_SYN);
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
		tau = TAUSYN;
	}

	for (i = 0; i < nn; i++) {
		
		if (i >= I_APP_NEURONS && USE_I_APP) {
			iapp = (time < I_APP_START || time > I_APP_END) ? I_APP : I_APP_STEP;	
		}
	
<<<<<<< HEAD
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
		
		if (self_connection) {
			if (arguments.divnn) {
				gsyn /= population;
			}
			current[I_S] =  gsyn * ((y[S] * (arguments.clustersize - 1))+ (y[P] * (population - arguments.clustersize))) * (y[V] - E_SYN);
		}

		else {
			current[I_S] = 0.0;
			if (arguments.delay >= arguments.stepsize) {
				for (j = 0; j < nn; ++j) {
					current[I_S] += weight[j + (i * nn)] * oldv[j]; //sums up products of weight and presynaptic y[S]
				}
			} 
			else {
				for (j = 0; j < nn; ++j) {
					current[I_S] += weight[j + (i * nn)] * y[S + (N * j)]; //sums up products of weight and presynaptic y[S]
				}
			}		
			//should the perturbation state variable thing be here instead?
			current[I_S] *= (arguments.divnn) ? ((gsyn / (nn - 1)) * (y[V + (N * i)] - E_SYN)) : (gsyn * (y[V + (N * i)] - E_SYN)); //multiplies synaptic current by maximum synaptic conductance and other stuff
			y[P + (N * i)] = current[I_S] ;	//sets perturbation state variable to synaptic current, doesn't affect simulation, purely for debugging purposes
		}
		
		//all of these (except for h) are using a method to prevent a divide by zero error I was encountering
		dydx[V + (N * i)] = (iapp - current[I_NA] - current[I_K] - current[I_M] - current[I_L] - current[I_S]) / CM;
		dydx[M + (N * i)] = ((0.32) * G((y[V + (N * i)] + 54), -4.) * (1. - y[M + (N * i)])) - ((0.28) * F((y[V + (N * i)] + 27), 5.) * y[M + (N * i)]);
		dydx[H + (N * i)] = 0.128 * exp(-(y[V + (N * i)] + 50.0) / 18.0) * (1.0 - y[H + (N * i)]) - 4.0 / (1.0 + exp(-(y[V + (N * i)] + 27.0) / 5.0)) * y[H + (N * i)];   
		dydx[NV + (N * i)] = ((0.032) * G((y[V + (N * i)] + 52), -5.) * (1. - y[NV + (N * i)])) - (0.5 * exp(-(y[V + (N * i)] + 57.) / 40.) * y[NV + (N * i)]); 
		dydx[MN + (N * i)] = ((3.209 * 0.0001) * G((y[V + (N * i)] + 30), -9.)  * (1.0 - y[MN + (N * i)])) + ((3.209 * 0.0001) * G((y[V + (N * i)] + 30), 9.) * y[MN + (N * i)]);
		
		if (self_connection) {
			if (arguments.delay >= arguments.stepsize) {
				dydx[S] = 2 * (1 + tanh(*oldv / 4.0)) * (1 - y[S]) - y[S] / tau; //uses *oldv which should be del, the delay pointer in the buffer
			}
			else {
				dydx[S] = 2 * (1 + tanh(y[V] / 4.0)) * (1 - y[S]) - y[S] / tau;
			}
			if (pertmode) {
				dydx[P] = 2 * (1 + tanh(*pert / 4.0)) * (1 - y[P]) - y[P] / tau;
			}
			else {
				dydx[P] = 2 * (1 + tanh((THRESHOLD) / 4.0)) * (1 - y[P]) - y[P] / tau;	//not sure why THRESHOLD is being used here?
			}
		}
		else {
			dydx[S + (N * i)] = ((2 * (1 + tanh(y[V + (i * N)] / 4.0))) * (1 - y[S + (i * N)])) - (y[S + (i * N)] / tau);
		}
	}
	return;
=======
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
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
}

void scan_(double *Y, int n, const char *filename) {
	FILE *fopen(),*sp;
	int i, j;
	sp = fopen(filename,"r");
	for (i = 0; i < n; ++i)
		if (fscanf(sp, "%lf\n", &Y[i]) != 1){
			fprintf(stderr, "Not enough variables in state.data file.\nNeed %d only %d given\n", N*nn, i-1);
			exit(1);
		}
	
	fclose(sp);
}

void dump_(double Y[]) {
	FILE *fopen(),*sp;
	int i, j, pos;
	sp = fopen("end.data","w");
	pos = 0;
	for (j = 0; j < nn; ++j) {
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
	//STU
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
		for (j = 0; j < (N * nn) - 1; j++) {
			fprintf(fp, "%f ", y[i][j]);
		}
		fprintf(fp, "%f\n", y[i][(N * nn) - 1]);
	}
	fclose(fp);
}

void makeallvolts(double** y, double *xx, int nstep, const char *filename) {
	int i, j;
	FILE *fopen(),*fp;
	fp = fopen(filename, "w");
	for (i = 0; i < nstep + 1; i++) {
		fprintf(fp, "%f ", xx[i]);
		for (j = 0; j < nn - 1; j++) {
			fprintf(fp, "%f ", y[i][V + (j * N)]);
		}
		fprintf(fp, "%f\n", y[i][V + (N * (nn - 1))]);
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

<<<<<<< HEAD
void snapshot(double** y, double *xx, int nstep, int var, double timestart, double timestop, Template *temp) {
	temp->steps = (int)(round((timestop - timestart) / arguments.stepsize));
=======

// This function is inefficient and should be optimized.
void snapshot(double** y, double *xx, int nstep, int var, double timestart, double timestop, Template *temp) {
	temp->steps = (int)(ceil((timestop - timestart) / STEPSIZE) + 2);
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
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
<<<<<<< HEAD
		temp->volts[temp->steps - 1] = STHRESHOLD;
=======
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
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
<<<<<<< HEAD
	fprintf(stderr, "This template contains %d steps.\n", temp->steps);
	fprintf(stderr, "Array of Voltages\n");
	printdarr(temp->volts, temp->steps);
	if (arguments.verbose) {
		fprintf(stderr, "Initial array of state variables\n");
		printdarr(temp->init, N);
		fprintf(stderr, "Initial array of buffer\n");
		printdarr(temp->ibuf, (int)(arguments.delay / arguments.stepsize));
	}
=======
	printf("This template contains %d steps.\n", temp->steps);
	printf("Array of Voltages\n");
	printdarr(temp->volts, temp->steps);
	printf("Initial array of state variables\n");
	printdarr(temp->init, N);
	printf("Initial array of buffer\n");
	printdarr(temp->ibuf, (int)(DELAY / STEPSIZE));
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
	
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

<<<<<<< HEAD
void pertsim(double normalperiod, Template spike, Phipair *trace, int tracedata, const char* tracename) {
	int i, k, nstep, targstep, flag, pertpos, psteps;
	double time, flag1, flag2;
	double v[N], vout[N], dv[N];			//v = variables (state and current), vout = output variables, dv = derivatives (fstate)
=======
int main() {
	//Variables to do with simulation
	long int nstep = (ENDTIME - STARTTIME) / STEPSIZE;	// This assumes the endtime-starttime is evenly divisible by the stepsize, which should always be true I think
	printf("The initial simulation will contain %d steps.\n", nstep);
	int i, k;
	double time;
	double v[N], vout[N], dv[N];					//v = variables (state and current), vout = output variables, dv = derivatives (fstate)
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
	double **y, *xx; 						//results variables, y[1..N][1..NSTEP+1], xx[1..NSTEP+1]
	extern double current[];				//external variable declaration
	
	//Variables to do with delay
<<<<<<< HEAD
	int dsteps = (int)(arguments.delay / arguments.stepsize);	//number of steps in the delay (i.e. number of elements in the buffer)
	double buf[dsteps];
	del = buf;
	int bufpos;
	
	psteps = (int)round(normalperiod / arguments.stepsize);
	nstep = (int)round((5.0 * normalperiod) / arguments.stepsize);
	
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
	copyab(spike.ibuf, buf, (int)(arguments.delay / arguments.stepsize));
	bufpos = spike.bufpos;
	del = &buf[bufpos]; //moves the pointer to the correct initial position in the buffer, unnecessary i believe
	pert = spike.volts;
	pertpos = 0;
	//~ printf("Template successfully initiated!\n");
	
	derivs(time, v, dv, del, &pweight);	//does running this one effect the phase/ perturbation? I don't think so but I'm not sure.
	rk4(v, dv, N, time, arguments.stepsize, vout, del, &pweight);
	
	xx[0] = time;
	for (i = 0; i < N; i++) {
		v[i] = vout[i]; 
		y[0][i] = v[i];
	}
	
	//~ printf("\nBeginning main for loop!\n\n");
	for (k = 0; k < nstep; k++) {			
		if (k == targstep && trace->phase >= 0.0) {	//activates perturbation mode if on correct step, allows derivs() to start using the "perturbation synapse" (a [pre-recorded stimulus of the same identical neuron)
				if (DBIT) {
					printf("k = %d\n", k);
				}
				pertmode = True;
				flag = True;
			}
		del = &buf[bufpos]; //moves the pointer one step ahead in the buffer
		derivs(time, v, dv, del, &pweight);
		rk4(v, dv, N, time, arguments.stepsize, vout, del, &pweight);
		*del = vout[0];
				
		time += arguments.stepsize;
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
					flag1 = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - arguments.stepsize))) + (time - arguments.stepsize);
					trace->fphi1 = ((double)(flag1) - (double)(normalperiod)) / (double)(normalperiod);
					//printf("The trace was done at phase %f. The step targeted was %d (time of flag1= %f), and the total number of steps in the unperturbed period was %d (%f ms).\n", trace->phase, targstep, time, psteps, normalperiod);
					//printf("f(phi)1 is %f.\n", trace->fphi1);
				}
				else {
					flag1 = k;
					trace->fphi1 = ((double)(k) - (double)(psteps))/ (double)(psteps);
				}
			}
		else if (flag && vout[0] >= THRESHOLD && v[0] < THRESHOLD && flag1 != 0.0 && flag2 == 0.0) {
				if (INTERPOLATE) {
					flag2 = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - arguments.stepsize))) + (time - arguments.stepsize);
					trace->fphi2 = ((double)(flag2) - (double)(flag1) - (double)(normalperiod))/ (double)(normalperiod);
				}
				else {
					flag2 = k;
					trace->fphi2 = ((double)(k) - (double)(flag1))/ (double)(psteps);
				}
			}
		for (i = 0; i < N; i++) {
			v[i] = vout[i];
			y[k + 1][i] = v[i];
		}
	}
	if (DBIT) {
 		printf("targstep = %d\n", targstep);
 	}
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

void template_init(Template *temp) {
	temp->ibuf = malloc(sizeof(double) * ((int)(arguments.delay / arguments.stepsize)));
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
	psteps = (int)round(normalperiod / arguments.stepsize);
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

void range_parser(int *start, int *end, const char *range) {
	char input[21];
	int i;
	char s[10];
	char e[10];
	strcpy(input, range);
	int len = strlen(range);
	for (i = 0; i < len; ++i) {
		if (input[i] == '-') {
			strcpy(e, &input[i + 1]);
			input[i] = '\0';
			strcpy(s, input);
			//printf("%s\n", s);
			//printf("%s\n", e);
			*start = atoi(s);
			*end = atoi(e);
			return;
		}
	}
	fprintf(stderr, "range could not be parsed. no '-' character was found\nenter ./lts --help for more information\n");
	exit(1);
	return;
}

void printargs(int argc, char **argv, const char *file) {
	int i;
	FILE *fopen(),*sp;
	sp = fopen(file, "w");
	for (i = 1; i < argc; ++i) {
		fprintf(sp, "%s ", argv[i]);
	}
	fprintf(sp, "\n");
	fclose(sp);
	return;
}

int main(int argc, char **argv) {
	

	/* Default values. */
	arguments.prc = 0;
	arguments.interval = 0;
	arguments.etime = 5000;
	arguments.lowprop = 0;
	arguments.lowstart = 0;
	arguments.lowend = 0;
	arguments.highprop = 0;
	arguments.highstart = 0;
	arguments.highend = 0;
	arguments.delay = 0;
	arguments.stepsize = 0.05;
	arguments.commit = 0;
	arguments.verbose = 0;
	arguments.graph = 0;
	arguments.divnn = 0;
	arguments.interval = 100;
	arguments.clustersize = 10;

	/* Parse our arguments; every option seen by parse_opt will
	be reflected in arguments. */
	argp_parse (&argp, argc, argv, 0, 0, &arguments);

	printf ("INPUTFILE = %s\nPRC = %s\nINTERVAL = %d\n", arguments.args[0], arguments.prc ? "yes" : "no", arguments.interval);
	printf("nn = %d\n", atoi(arguments.args[1]));
	if (arguments.divnn) {
		printf("DIVNN IS enabled\n");
	}
	else {
		printf("DIVNN IS NOT enabled\n");
	}

	nn = atoi(arguments.args[1]);
	population = nn;
	iapps = malloc(sizeof(double) * nn);
	
	gsyn = G_SYN;
	tau = TAUSYN;
	//Variables to do with simulation
	long int nstep = (arguments.etime - STARTTIME) / arguments.stepsize;	// This assumes the endtime-starttime is evenly divisible by the stepsize, which should always be true I think
	fprintf(stderr, "The initial simulation will contain %d steps.\n", nstep);
	int i, j, k, pre;
	double time;
	double v[N * nn], vout[N * nn], dv[N * nn];					//v = variables (state and current), vout = output variables, dv = derivatives (fstate)
	double **y, *xx; 						//results variables, y[1..N][1..NSTEP+1], xx[1..NSTEP+1]
	extern double current[];				//external variable declaration
	double weight[nn*nn];
	scan_(weight, nn*nn, "weights.data");
	
	fprintf(stderr, "DELAY = %fms.\n", arguments.delay);
	if (arguments.lowprop) {
		fprintf(stderr, "Using Low Dose Propofol starting at time %d and ending at %d.\n", arguments.lowstart, arguments.lowend);
	}
	if (arguments.highprop) {
		fprintf(stderr, "Using High Dose Propofol starting at time %d and ending at %d.\n", arguments.highstart, arguments.highend);
	}
	/*
	fprintf(stderr, "Printing data to file ");
	fprintf(stderr, FULLNAME);
	fprintf(stderr, "\n");
	*/

	//Variables to do with delay
	int dsteps = (int)(arguments.delay / arguments.stepsize);	//number of steps in the delay (i.e. number of elements in the buffer)
	double* buf[dsteps];
	for (i = 0; i < dsteps; ++i) {
		buf[i] = (double*) malloc(sizeof(double) * nn);		//allocates memory for an array of arrays, length depending on the number of neurons in the simulation
	}
	int bufpos = 0;								//holds the position in the buffer that the del  pointer is at
	del = buf[bufpos];
	

=======
	int dsteps = (int)(DELAY / STEPSIZE);	//number of steps in the delay (i.e. number of elements in the buffer)
	double buf[dsteps];
	del = buf;
	int bufpos;								//holds the position in the buffer that the del  pointer is at
	
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
	//Variables to do with initial PRC measurements
	double normalperiod;					//unperturbed period of the oscillation, difference between snd_time and fst_time;
	int psteps;								//number of steps in the unperturbed period
	double sptimes[SAMPLESIZE + OFFSET];	//array of times of sp(PRCSKIP) ? ((((double)(k) - (double)(psteps)) - (double)(psteps)) / (double)(psteps)) : ((double)(k) - (double)(psteps))/ (double)(psteps);iking, differences will be averaged to find unperturbed period of oscillation
	int spikecount = 0;						//holds the location the simulation has reached in sptimes
	double spdiffs[SAMPLESIZE - 1];			//array of differences in the times of spiking, averaged to find the normalperiod
	double sumdiffs = 0;					//holds the sum of differences in times of sptimes, used for averaging
	double periodsd;						//standard deviation of the averaged periods
	
<<<<<<< HEAD
	//Allocating memory for the storage arrays, checking if I can, so that I don't run out of memory
	y = (double**) malloc(sizeof(double*) * (nstep + 1));
	for (i = 0; i < (nstep + 1); i++) {
		y[i] = (double*) malloc(sizeof(double) * (N * nn));
=======
	//Variables for storing spike snapshot
	double fthresh = -1.0;					//time at which voltage crosses threshold on the way up
	double sndthresh = -1.0;				//time at which voltage crosses threshold on the way down
	Template spike;

	
	//Allocating memory for the storage arrays, checking if I can, so that I don't run out of memory
	y = (double**) malloc(sizeof(double*) * (nstep + 1));
	for (i = 0; i < (nstep + 1); i++) {
		y[i] = (double*) malloc(sizeof(double) * N);
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
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
<<<<<<< HEAD
	
	time = STARTTIME;
	scan_(v, N*nn, "state.data");				//scanning in initial variables (state variables only) 

	if (arguments.delay >= arguments.stepsize) {
		for (i = 0; i < (dsteps); ++i) {//sets every double in buffer(s) to be equal to 0
			for (j = 0; j < nn; ++j) {
					buf[i][j] = 0.0;
			}
		}
	}

=======
	
	//Variables to do with conducting prc
	extern int pertmode;
	extern double *pert;
	
	
	
	time = STARTTIME;
	scan_(v);				//scanning in initial variables (state variables only) 
	
	for (i = 0; i < (dsteps); ++i) {//sets every double in buffer to be equal to the steady state (initial) voltage that was just scanned in
		buf[i] = v[0];
	}
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
	
	iapps[0] = I_APP;
	for (i = 1; i < nn; ++i) {						//adding heterogeneity through slightly different iapps like mccarthy does
		iapps[i] = iapps[i - 1] + 0.005;
	}
	
	derivs(time, v, dv, del, weight);
	rk4(v, dv, (N * nn), time, arguments.stepsize, vout, del, weight);
	
	
	xx[0] = STARTTIME;		
	for (i = 0; i < (N * nn); ++i) {
		v[i] = vout[i]; 
		y[0][i] = v[i];
		
	}

	for (k = 0; k < nstep; ++k) {
		
		if (arguments.lowprop && time > arguments.lowstart && time < arguments.lowend) {	//changes gsyn to match the correct level of propofol in the simulation
			gsyn = (time > arguments.lowstart && time < arguments.lowend) ? LOWPROP_GSYN : G_SYN; 
			tau  = (time > arguments.lowstart && time < arguments.lowend) ? LOWPROP_TAU : TAUSYN; 
		}
		else if (arguments.highprop && time > arguments.highstart && time < arguments.highend) {
			gsyn = (time > arguments.highstart && time < arguments.highend) ? HIGHPROP_GSYN : G_SYN;
			tau  = (time > arguments.highstart && time < arguments.highend) ? HIGHPROP_TAU : TAUSYN;
		}
		else {	
			gsyn = (G_SYN);
			tau  = TAUSYN;
		}
		
		del = buf[bufpos]; //moves the pointer one step ahead in the buffer
		derivs(time, v, dv, del, weight);										//actually does the step
		rk4(v, dv, (N * nn), time, arguments.stepsize, vout, del, weight);				//actually does the step
		
		if (arguments.delay >= arguments.stepsize) {
			for (j = 0; j < nn; ++j) {								
				del[j] = vout[S + (N * j)];							//dereferences the delay pointer, and puts the previous y[S] in the same spot
			}
		}
		if (vout[0] >= THRESHOLD && v[0] < THRESHOLD) {
			if (spikecount < (SAMPLESIZE + OFFSET)) {
				if (INTERPOLATE) {
<<<<<<< HEAD
					sptimes[spikecount] = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - arguments.stepsize))) + (time - arguments.stepsize);
=======
					sptimes[spikecount] = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - STEPSIZE))) + (time - STEPSIZE);
					printf("sptimes[spikecount] with interpolation is %f. Without is %f. The difference is %f.\n", sptimes[spikecount], time, (sptimes[spikecount] - time));
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
				}
				else {
					sptimes[spikecount] = time;
				}
			}
			++spikecount;			//incremented at the end so it can be used as position in sptimes			
		}
		
<<<<<<< HEAD
		time += arguments.stepsize;
=======
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
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
		xx[k + 1] = time;
		
		if (arguments.delay >= arguments.stepsize) {
			bufpos = (++bufpos)%dsteps;
		}
		
		for (i = 0; i < (N * nn); i++) {
			v[i] = vout[i];
			y[k + 1][i] = v[i];
		}
	}
<<<<<<< HEAD
	
	makeallvolts(y, xx, nstep, FULLNAME);
	
	if (arguments.verbose) {
		dump_(vout);
		makedata(y, xx, nstep, V, "v.data");
		//makedata(y, xx, nstep, (V + N), "v2.data");
		makedata(y, xx, nstep, M, "m.data");
		makedata(y, xx, nstep, H, "h.data");
		makedata(y, xx, nstep, NV, "n.data");
		makedata(y, xx, nstep, S, "s.data");
		makefull(y, xx, nstep, "full.data");
		
		
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
		//makefullsingle(y, xx, nstep, 18, "18full.data");
		//~ makefullsingle(y, xx, nstep, 19, "19full.data");
	}
	else {
		fprintf(stderr, "\n\nSince verbose == 0, v-n.data are not being written Only vfull.data\n\n");
	}
	fprintf(stderr,"This simulation counted %d spikes of Neuron[0].\n", spikecount);
=======
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
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
	if (spikecount >= (SAMPLESIZE + OFFSET)) {
		for (i = OFFSET; i < SAMPLESIZE + OFFSET - 1; ++i) {		//calculates differences between spike times to find each period
			sumdiffs += sptimes[i + 1] - sptimes[i];
			spdiffs[i - OFFSET] = sptimes[i + 1] - sptimes[i];
		}
<<<<<<< HEAD
		if (arguments.verbose) {
			printperiod(spdiffs, SAMPLESIZE - 1, "period.data");
		}
		normalperiod = sumdiffs / (SAMPLESIZE - 1);
		psteps = (int)round(normalperiod / arguments.stepsize);
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
	free(iapps);

	//Starting PRC simulation, if requested
	if (arguments.prc || DO_TRACE) {

		nn = 1;
		self_connection = 1;
		printf("Running PRC on clustersize %d!\nGraph the main simulation to see whether that's appropriate for the PRC being run.\n", arguments.clustersize);
		
		//Variables to do with delay
		int dsteps = (int)(arguments.delay / arguments.stepsize);	//number of steps in the delay (i.e. number of elements in the buffer)
		double prcbuf[dsteps];
		del = prcbuf;
		int bufpos = 0;

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
		template_init(&spike);

		int startstep;
		fprintf(stderr, "\n\n\n\nDoing PRC and/or Trace\n-running preliminary PRC simulation\n");

		//Allocating memory for the storage arrays, checking if I can, so that I don't run out of memory
		y = (double**) malloc(sizeof(double*) * (nstep + 1));
		for (i = 0; i < (nstep + 1); i++) {
			y[i] = (double*) malloc(sizeof(double) * N);
			if (y[i] == NULL) {
				fprintf(stderr, "Ran out of memory for storage array at y[%d]\n", i);
				return 0;
			}
		}
		xx = (double*) malloc(sizeof(double) * (nstep + 1));
		if (xx == NULL) {
			fprintf(stderr, "Ran out of memory for storage array xx]\n");
			return 0;
		} 
		
		//Variables to do with conducting prc
		extern int pertmode;
		extern double *pert;
	
		time = STARTTIME;
		scan_(v, N, "state.data");				//scanning in initial variables (state variables only)

		fprintf(stderr, "Here are the starting variables for the single self-connected prc neuron.\n");
		printdarr(v, N);

		fprintf(stderr, "-scanning in first 7 variables of state.data\n");
		fprintf(stderr, "-this PRC might not be accurate unless the initial conditions in the first 7 are representative of all the others\n"); 
		
		for (i = 0; i < (dsteps); ++i) {//sets every double in buffer to be equal to the steady state (initial) voltage that was just scanned in
			prcbuf[i] = v[0];
		}
		
		derivs(time, v, dv, del, &pweight);
		
		rk4(v, dv, N, time, arguments.stepsize, vout, del, &pweight);
		
		xx[0] = STARTTIME;		
		for (i = 0; i < N; i++) {
			v[i] = vout[i]; 
			y[0][i] = v[i];
		}

		for (k = 0; k < nstep; k++) {
			del = &prcbuf[bufpos]; //moves the pointer one step ahead in the buffer
			derivs(time, v, dv, del, &pweight);
			rk4(v, dv, N, time, arguments.stepsize, vout, del, &pweight);
			*del = vout[0];
			
			if (vout[0] >= THRESHOLD && v[0] < THRESHOLD) {
				if (spikecount < (SAMPLESIZE + OFFSET)) {
					if (INTERPOLATE) {
						sptimes[spikecount] = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - arguments.stepsize))) + (time - arguments.stepsize);
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
						fthresh = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - arguments.stepsize))) + (time - arguments.stepsize);
					}
					else {
						fthresh = time;
					}
					spike.init[0] = -50.0;
					printf("Using variables at time %f (interpolated time %f) for start of spike template.\n", time, fthresh);
					if (DBIT) {
						printf("Current state variables, except voltage is interpolated to 0.\n");
						printdarr(vout, N);
					}
					startstep = k + 1;
					
					for (i = 1; i < N; ++i) {			//puts initial state variables into spike template
						spike.init[i] = vout[i];
					}
					if (DBIT) {
						printf("Printing the spike.init array just in case.\n");
						printdarr(spike.init, N);
					}
					
					printf("About to copy buffer and buffer position to template.\n");
					if (DBIT) {
 						printf("Buffer array.\n");
						printdarr(prcbuf, dsteps);
						printf("bufpos = %f\n", bufpos);
					}
					
					for (i = 0; i < dsteps; ++i) {
						spike.ibuf[i] = prcbuf[i];			//puts initial buffer into spike template
					}
					spike.bufpos = bufpos;
					if (DBIT) {
						printf("Done copying. spike.ibuf\n");
						printdarr(spike.ibuf, dsteps);
						printf("spike.bufpos = %f\n", spike.bufpos);
					}
				}
				else if (fthresh != -1.0 && sndthresh == -1.0 && vout[0] <= STHRESHOLD && v[0] > STHRESHOLD) {
					if (INTERPOLATE) {
						sndthresh = ((((THRESHOLD) - v[0]) / (vout[0] - v[0])) * (time - (time - arguments.stepsize))) + (time - arguments.stepsize);
					}
					else {
						sndthresh = time;
					}
				}
			}
					
			time += arguments.stepsize;
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
		if (arguments.verbose) {
			makedata(y, xx, nstep, V, "v.data");
			makedata(y, xx, nstep, M, "m.data");
			makedata(y, xx, nstep, H, "h.data");
			makedata(y, xx, nstep, NV, "n.data");
		}
		else {
			fprintf(stderr, "\n\nSince PLONG == 0, v-n.data are not being written\n\n");
		}
		if (arguments.verbose) {
			dump_(vout);
		}
		fprintf(stderr, "This simulation counted %d spikes in all.\n", spikecount);
		if (spikecount >= (SAMPLESIZE + OFFSET)) {
			for (i = OFFSET; i < SAMPLESIZE + OFFSET - 1; ++i) {		//calculates differences between spike times to find each period
				sumdiffs += sptimes[i + 1] - sptimes[i];
				spdiffs[i - OFFSET] = sptimes[i + 1] - sptimes[i];
			}
			//printperiod(spdiffs, SAMPLESIZE - 1, "period.data");
			normalperiod = sumdiffs / (SAMPLESIZE - 1);
			psteps = (int)round(normalperiod / arguments.stepsize);
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
		fprintf(stderr, "fthresh = %f and sndthresh = %f\n", fthresh, sndthresh);
		
		snapshot(y, xx, nstep, V, fthresh, sndthresh, &spike);
		printemp(&spike);
		
		if (arguments.prc) {
			prc(spike, arguments.interval, normalperiod);
		}
		if (DO_TRACE) {
			Phipair test;
			test.phase = -1.0;
			if (DBIT) {
				fprintf(stderr, "test.fphi1 = %f\n", test.fphi1);
			}
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
			if (DBIT) {
				fprintf(stderr, "test.fphi1 = %f\n", test.fphi1);	
				fprintf(stderr, "trace.fphi1 = %f\n", trace.fphi1);
			}	
		}
		for (i = 0; i < (nstep + 1); i++) {		
			free(y[i]);
		}
		free(xx);	
		free(spike.volts);
		free(spike.ibuf);
		

		if (arguments.graph) {
			char test[300];
			sprintf(test, "python2 testline.py prc1.data %f %f", arguments.delay, normalperiod);
			fprintf(stderr, "%s\n", test);
			system(test);
		}
	}	
=======
		printdarr(spdiffs, SAMPLESIZE - 1);
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
		printf("v.data-n.data as well as end.data were still written, but no trace or prc processes occured.");
		return 0;
	}
	printf("fthresh = %f and sndthresh = %f\n", fthresh, sndthresh);
	
	snapshot(y, xx, nstep, V, fthresh, sndthresh, &spike);
	printemp(&spike);
	
	int prcsteps;
	extern int prcmode;
	double phase;
	prcsteps = psteps * 5;
	double targphase;
	int targstep;
	int pertpos;
	int flag;		//perturbation happened
	int flag1;		//first period complete after perturbation
	int flag2;		//second period complete after perturbation
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
			if (flag && vout[0] >= THRESHOLD && v[0] < THRESHOLD && !flag1) {
				flag1 = k;
				trace.fphi1 = (PRCSKIP) ? ((((double)(k) - (double)(psteps)) - (double)(psteps)) / (double)(psteps)) : ((double)(k) - (double)(psteps))/ (double)(psteps);
				printf("The trace was done at phase %f. The step targeted was %d (time = %f), and the total number of steps in the unperturbed period was %d.\n", TPHASE, targstep, time, psteps);
				printf("f(phi)1 is %f.\n", trace.fphi1);
			}
			else if (flag && vout[0] >= THRESHOLD && v[0] < THRESHOLD && flag1 != 0 && !flag2) {
				flag2 = k;
				trace.fphi2 = ((double)(k) - (double)(flag1))/ (double)(psteps);
				//~ printf("k %d ptime %f targstep %d prc[prcpos].phase %f time %f prc[prcpos].fphi1 %f flag1 %d flag2 %d prc[prcpos].fphi2 %f\n", k, ptime, targstep, prc[prcpos].phase, time, prc[prcpos].fphi1, flag1, flag2, prc[prcpos].fphi2);
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
			flag1 = False;
			flag2 = False;
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
				if (flag && vout[0] >= THRESHOLD && v[0] < THRESHOLD && flag1 == 0) {
					flag1 = k;
					prc[prcpos].fphi1 = (PRCSKIP) ? ((((double)(k) - (double)(psteps)) - (double)(psteps)) / (double)(psteps)) : ((double)(k) - (double)(psteps))/ (double)(psteps);
					prc[prcpos].phase = targphase;
					//~ ((double)targstep - (double)psteps) / (k - psteps);
					
				}
				else if (flag && vout[0] >= THRESHOLD && v[0] < THRESHOLD && flag1 != 0 && !flag2) {
					flag2 = k;
					prc[prcpos].fphi2 = ((double)(k) - (double)(flag1))/ (double)(psteps);
					printf("k %d ptime %f targstep %d prc[prcpos].phase %f time %f prc[prcpos].fphi1 %f flag1 %d flag2 %d prc[prcpos].fphi2 %f\n", k, ptime, targstep, prc[prcpos].phase, time, prc[prcpos].fphi1, flag1, flag2, prc[prcpos].fphi2);
					prcpos++;
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
	
	for (i = 0; i < (nstep + 1); i++) {		
		free(y[i]);
	}
		
	
	free(spike.volts);
>>>>>>> 9e28c4720c333de062b7131c3b4d59240c94fce0
	
	printargs(argc, argv, "args.txt");
	if (arguments.commit) {
		system("git commit -a");
	}
	return 0;
}
