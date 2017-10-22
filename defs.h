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
#define STEPSIZE 0.05
#define DELAY 3.5			//delay must evenly divide stepsize, and it is only used if it is >= stepsize
#define THRESHOLD -50.0		//the voltage at which it counts a spike has occured, used to measure both nonperturbed and perturbed period for PRC
#define STHRESHOLD -50.0	//threshold used to measure just the spike, not the period between spikes
#define SAMPLESIZE 5 		//number of spikes that are averaged together to give unperturbed period
#define OFFSET 10			//number of spikes that are skipped to allow the simulation to "cool down" before it starts measuring the period
#define POPULATION 20		//number of neurons in the whole population
#define MYCLUSTER 10		//number of neurons in the simulated neuron's population
#define True 1
#define False 0
#define INTERPOLATE 0
#define PLONG 1
#define FULLNAME "lowhigh.data"
#define DBIT 1
#define DIVNN 1
#define G(X,Y) ( (fabs((X)/(Y))<1e-6) ? ((Y)*((X)/(Y)/2. - 1.)) : ((X)/(1. - exp( (X)/ (Y) ))) )
#define F(X,Y) ( (fabs((X)/(Y))<1e-6) ? ((Y)*(1.-(X)/(Y)/2.)) : ((X)/(exp( (X)/ (Y) ) -1)) )
#define PERTENDTIME 5000	//separate endtime for prc stuff in order to differentiate it from main simulation
