//prcfuncs
double current[C];	//external current variable, similar to how Canavier did it
static double *del;
double gsyn, tau;
double iapps[NN];
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
	double iapp;
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
		
		current[I_S] *= (DIVNN) ? ((gsyn / (NN - 1)) * (y[V + (N * i)] - E_SYN)) : (gsyn * (y[V + (N * i)] - E_SYN)); //multiplies synaptic current by maximum synaptic conductance and other stuff
		
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

void prcderivs(double time, double *y, double *dydx, double *oldv) { 
	double iapp, gsyn, tau;
	extern double *pert;
	
	
	if (USE_I_APP && !(prcmode)) {
		iapp = (time < I_APP_START || time > I_APP_END) ? I_APP : I_APP_STEP;
		
	}
	else {
		iapp = I_APP;
	}
	
	if (USE_LOWPROPOFOL) {
		gsyn = ((time > LOW_PROPOFOL_START && time < LOW_PROPOFOL_END) || prcmode) ? LOWPROP_GSYN : G_SYN;
		tau = ((time > LOW_PROPOFOL_START && time < LOW_PROPOFOL_END) || prcmode) ? LOWPROP_TAU : TAUSYN; 
	}
	else if (USE_HIGHPROPOFOL) {
		gsyn = ((time > HIGH_PROPOFOL_START && time < HIGH_PROPOFOL_END) || prcmode) ? HIGHPROP_GSYN : G_SYN;
		tau = ((time > HIGH_PROPOFOL_START && time < HIGH_PROPOFOL_END) || prcmode) ? HIGHPROP_TAU : TAUSYN;
	}
	else {	
		gsyn = G_SYN;
		tau = TAUSYN;
	}
	
	if (DIVNN) {
		gsyn /= POPULATION;
	}
	
	current[I_NA] = G_NA * y[H] * y[M] * y[M] * y[M] * (y[V] - E_NA);
	current[I_K] =  G_K * y[NV] * y[NV] * y[NV] * y[NV] * (y[V] - E_K);
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

void prcrk4(double y[], double dydx[], int n, double x, double h, double yout[], double *oldv) {
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
	
	prcderivs(xh, yt, dyt, oldv);				//second step
	
	for (i = 0; i < n; i++) {
		yt[i] = y[i] + hh * dyt[i];
	}
	
	prcderivs(xh, yt, dym, oldv);				//third step
	
	for (i = 0; i < n; i++) {
		yt[i] = y[i] + h * dym[i];
		dym[i] += dyt[i];
	}
	
	prcderivs(x + h, yt, dyt, oldv);			//fourth step
	
	for (i = 0; i < n; i++) {			//Accumulate increments with proper weights.
		yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
	}
	
	free(dym);
	free(dyt);
	free(yt);
}

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
	
	prcderivs(time, v, dv, del);	//does running this one effect the phase/ perturbation? I don't think so but I'm not sure.
	prcrk4(v, dv, N, time, STEPSIZE, vout, del);
	
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
		prcderivs(time, v, dv, del);
		prcrk4(v, dv, N, time, STEPSIZE, vout, del);
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



















