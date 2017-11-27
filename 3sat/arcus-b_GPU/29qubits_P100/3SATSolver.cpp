/** @file 
 * Simulates adiabatic 3SAT solver
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <sys/types.h>   		// needed on ARC to find ssize_t
#include <string.h>

#include "3SATGenerator.h"		// new SATs aren't generated here though
# include "../../../QuEST_GPU_v0.5.0/precision.h"
# include "../../../QuEST_GPU_v0.5.0/qubits.h"
# include <sys/time.h>

#define CLAUSE_SIZE 3
#define CLAUSE_BOOL_RATIO 4.267

#define N_TRIALS 3

const char *OUTPUT_FORMAT = "%.5e";

const long double PI		=3.14159265358979323846264338327950288419716939937510;
const long double SQRT2	= 1.41421356237309504880168872420969807856967187537694;

QuESTEnv env;

int numLocal=0, numDist=0;

REAL system_timer (void) {
        struct timeval time;
        gettimeofday (&time, NULL);
        return time.tv_sec + time.tv_usec / 1000000.0;

}

void saveInteger(int value, char *filename, char *varname) {
	
	FILE *file = fopen(filename, "a");
	fprintf(file, "%s=%d;\n\n", varname, value);
	fclose(file);
}


void saveDouble(double value, char *filename, char *varname) {
	
	FILE *file = fopen(filename, "a");
	fprintf(file, "%s=%.10f;\n\n", varname, value);
	fclose(file);
}


void saveIntegerArray(int *array, int length, char *filename, char *varname) {
	
	FILE *file = fopen(filename, "a");
	
	// open mathematica brace
	fprintf(file, varname);
	fprintf(file, "={");
	for (int i=0; i < length; i++) {
		fprintf(file, "%d", array[i]);
		if (i < length-1)
			fprintf(file, ",");
	}
	
	// close mathematica brace (and file)
	fprintf(file, "};\n\n");
	fclose(file);
}


void saveDoubleArray(double *array, int length, char *filename, char *varname) {
	
	FILE *file = fopen(filename, "a");
	
	// open mathematica brace
	fprintf(file, varname);
	fprintf(file, "={");
	
	for (int i=0; i < length; i++) {
		//fprintf(file, "%.5e,", array[i]);
		
		// get 'baseE' format string in buffer
		ssize_t numSize = snprintf(NULL, 0, OUTPUT_FORMAT, array[i]);
		int numChars = numSize/sizeof(char);
		char *numBuffer = (char*) malloc(numSize + 1);
		snprintf(numBuffer, numSize+1, OUTPUT_FORMAT, array[i]);
		
		// rewrite 'baseE' as 'base10'
		for (int j=0; j < numChars; j++)
			if (numBuffer[j] == 'e')
				fprintf(file, "*10^");
			else
				fprintf(file, "%c", numBuffer[j]);
		free(numBuffer);
		
		// separate elements with commas
		if (i < length-1)
			fprintf(file, ",");
		
	}
	
	// close mathematica brace (and file)
	fprintf(file, "};\n\n");
	fclose(file);
}


int recordProgress(char *filename, int cycle, int numCycles, int lastCycleTime) {
	
	FILE *file = fopen(filename, "a");
	time_t timeNow = time(0);
	if (cycle > 0) 
		fprintf(file, "took %ds\n", (int) timeNow - lastCycleTime);
	
	fprintf(file, "%d/%d ", cycle+1, numCycles);
	fclose(file);
	return timeNow;
}


void createOrClearFile(char *filename) {
	FILE *file = fopen(filename,"w");
	fclose(file);
}


void extractSolution(MultiQubit qubits, int *solBits, REAL *solProb) {
	int solInt = 0;
	getLargestProbEl(qubits, solProb, &solInt);
	convertToBinary(solInt, qubits.numQubits, solBits);
}


void transformQubit(MultiQubit qubits, int qubitInd, double theta, int qubitBool, double *angles) {
	
	double newAngle = (qubitBool)? theta : -theta;
	double oldAngle = angles[qubitInd];
	
	// if angles angree, no transformation required
	if (oldAngle == newAngle)
		return;
		
	// update angle record
	angles[qubitInd] = newAngle;
	
	// transform qubit from old to new angle
	Complex alpha, beta;
	alpha.imag = beta.imag = 0;
	alpha.real = cos(.5*(oldAngle - newAngle));
	beta. real = sin(.5*(oldAngle - newAngle));
	rotateQubit(qubits, qubitInd, alpha, beta);

	if (qubitInd==28) numDist++;
	else numLocal++;
}


int *solve3SAT(

		// inputs
		int* equ, 
		int numBools, int numCycles, int numClauses, double *thetaFactors,
		
		// outputs
		double *totalSolProb, int *runTime, double **qubitZeroProbs,
		double **clausePassProbs, double **clauseReachedFailProbs,
		
		// progress
		int recordProgressFlag, char *progressFN
	) {

	/*
	 * PREPARE SIMULATOR
	 */
	 
	// load QuEST, allocate qubits
	MultiQubit qubits; 
	createMultiQubit(&qubits, numBools, env);

	reportQuESTEnv(env);
	reportMultiQubitParams(qubits);
	
	// prepare data collection structures
	double solProb;				// probability of measuring sol in final state
	double totalProb;			// total probability of reaching final state
	//*totalSolProb;				// prob of reaching final state and measuring sol
	//*qubitZeroProbs;			// prob of each qubit being zero after clause cycle
	//*clausePassProbs;			// probs of passing clauses, given reached
	//*clauseReachedFailProbs;	// probs of reaching then failing clauses
	//*runTime;					// algorithm execution time (in seconds)
	
	// prepare iteration vars
	int cycleNum;		// number of clause-check-cycles so far
	int clauseInd;		// index of the clause being checked
	int localInd;		// iterates 0 to CLAUSE_SIZE
	int termInd;		// index of the considered term in equ
	int termQb;			// index of the considered term's qubit in register
	int termBool;		// 0 if considered term appears negated, else 1
	int clauseProbInd;	// index of the considered clause in prob arrays
	double probCheck; 	// probability of successful clause check
	double thetaFactor;	// the current factor theta/(PI/2); adiabatic variable
	double theta;		// the adiabatic variable theta, controlling measurement
	double *qubitAngles;	// the 'angle' qubit is currently rotated to
	
	// allocate needed memory
	totalProb = 1;
	int *sol				= (int*) malloc(numBools * sizeof(int));
	*qubitZeroProbs			= (double*) malloc(numBools   * numCycles * sizeof(double));
	*clausePassProbs			= (double*) malloc(numClauses * numCycles * sizeof(double));
	*clauseReachedFailProbs	= (double*) malloc(numClauses * numCycles * sizeof(double));
	qubitAngles 			= (double*) malloc(numBools * sizeof(double));
	
	// initial angles are -pi/2, corresponding to I (no previous gates)
	for (int i=0; i < numBools; i++) {
		qubitAngles[i] = -PI/2;
	}
	
	// prepare progress reporting
	time_t lastCycleTime = time(0);
	if (recordProgressFlag)
		createOrClearFile(progressFN);

	/*
	 * PERFORM ALGORITHM
	 */
	 
	// start timeing execution
	time_t startTime = time(0);
	 
	// populate all states
	initStatePlus(&qubits);
	
	// repeat many clause-check cycles
	for (cycleNum=0; cycleNum<numCycles; cycleNum++) {
		// optionally record progress
		if (recordProgressFlag)
			lastCycleTime = recordProgress(
				progressFN, cycleNum, numCycles, lastCycleTime);
		
		// update theta
		thetaFactor = thetaFactors[cycleNum];
		theta = thetaFactor * PI/2;
		if (env.rank==0) printf("cycle %d/%d, theta=%f pi/2\n", cycleNum+1, numCycles, thetaFactor);
		
		// check each clause
		for (clauseInd=0; clauseInd<numClauses; clauseInd++) {
			// transform each qubit-in-clause
			for (localInd=0; localInd<CLAUSE_SIZE; localInd++) {
				
				// locate qubit
				termInd = clauseInd*CLAUSE_SIZE + localInd;
				termQb = abs(equ[termInd]) - 1;
				termBool = (equ[termInd] > 0);
				
				// transform into previous angle, then out of undesired to |0>
				transformQubit(qubits, termQb, theta, termBool, qubitAngles);
			}
			
			// remove undesired state
			probCheck = filterOut111(
				qubits, 
				abs(equ[termInd-2]) - 1,
				abs(equ[termInd-1]) - 1,
				abs(equ[termInd-0]) - 1);
				
			// process prob statistics
			clauseProbInd = cycleNum*numClauses + clauseInd;
			(*clausePassProbs)[clauseProbInd] = probCheck;
			(*clauseReachedFailProbs)[clauseProbInd] = totalProb*(1-probCheck);
			totalProb *= probCheck;
		}
		// collect zero projection of each qubit
		for (int qbInd=0; qbInd < numBools; qbInd++) {
			int qbClauseInd = numBools*cycleNum + qbInd;
			double zeroProb = findProbabilityOfOutcome(qubits, qbInd, 0);
			(*qubitZeroProbs)[qbClauseInd] = zeroProb;
		}
	}

	// rotate all qubits to the final -pi/2 angle (corresponding to I)
	for (int qbInd=0; qbInd < numBools; qbInd++)
		transformQubit(qubits, qbInd, -PI/2, 1, qubitAngles);
	
	/*
	 * PROCESS RESULTS
	 */
	
	// extract solution state from final state
	extractSolution(qubits, sol, &solProb);
	*totalSolProb = solProb * totalProb;
	

	// stop timing execution
	*runTime = (int)(time(0) - startTime);
	
	/*
	 * FREE MEMORY
	 */
	
	destroyMultiQubit(qubits, env); 
	
	free(qubitAngles);

	return sol;
}


int main (int narg, char** varg) {

	initQuESTEnv(&env);

	/*
	 * GET PARAMETERS
	 */
	
	// 3SATSolver inFN, outFN numCycles initThetaFactor updateTheta[ progressFN]
	if (narg != 6 && narg != 7) {
		printf("ERROR: call as ./3SATSolver "
			   "inFN outFN numCycles "
			   "initThetaFactor updateTheta[ progressFN]\n");
		closeQuESTEnv(env);
		return EXIT_FAILURE;
	}
	
	char *inFN = (char*) malloc(strlen(varg[1])+1);
	char *outFN = (char*) malloc(strlen(varg[2])+1);
	strcpy(inFN, varg[1]);
	strcpy(outFN, varg[2]);
	
	int numCycles = atoi(varg[3]);
	double initThetaFactor;
	sscanf(varg[4], "%lf", &initThetaFactor);
	int updateTheta = atoi(varg[5]);
	
	char *progressFN;
	int progressFlag = (narg == 7);
	if (progressFlag) {
		progressFN = (char*) malloc(strlen(varg[6])+1);
		strcpy(progressFN, varg[6]);
	} else
		progressFN = (char*) malloc(sizeof(char));
	
	
	/*
	 * PREPARE NEW 3SAT PROBLEM
	 */
	
	/*
	 // generate random problem
	int numBools = 20;
	int numClauses = getNumClauses(numBools);
	int *equ = malloc(numClauses * CLAUSE_SIZE * sizeof(int));
	int *sol = malloc(numBools * sizeof(int));
	printf("generating %dqb 3SAT problem\n", numBools);
	getRandomEquAndSol(equ, sol, numBools, numClauses);
	saveEquation(filename, equ, numClauses);
	*/
	
	
	/*
	 * LOAD EXISTING 3SAT PROBLEM
	 */
	 
	int numBools;
	int numClauses;
	int *equ = loadEquation(inFN, &numBools, &numClauses);

	
	/*
	 * SOLVE 3SAT PROBLEM
	 */
	 
	// collect theta values
	double *thetaFactors = (double*) malloc(numCycles * sizeof(double));
	for (int cycleNum=0; cycleNum < numCycles; cycleNum++) {

		if (updateTheta) {
			long double progress = cycleNum / (double) numCycles;
			double thetaFactor = initThetaFactor + (1-initThetaFactor)*pow(progress, 3);
			thetaFactors[cycleNum] = thetaFactor;
		} else
			thetaFactors[cycleNum] = initThetaFactor;
	}
	
	// declare outputs
	int runTime;
	double totalSolProb;
	double *qubitZeroProbs;
	double *clausePassProbs;
	double *clauseReachedFailProbs;

	// prepare timing
	REAL wtime_start, wtime_stop;
	REAL *timingVec;
	FILE *timing;
	int trial;
	char envString[255];
	char filename[255];
	if (env.rank==0) timingVec = (REAL*) malloc(N_TRIALS*sizeof(timingVec));

	if (env.rank==0){
		MultiQubit testMQ;
		createMultiQubit(&testMQ, numBools, env);
		getEnvironmentString(env, testMQ, envString);
		destroyMultiQubit(testMQ, env);
		sprintf(filename, "TIMING3SAT_%s.csv", envString);
		timing = fopen(filename, "w");
		fprintf(timing, "time(s), standardDev, maxDelta, minDelta\n");
	}


	int *sol;
	for (trial=0; trial<N_TRIALS; trial++){
		// start timing	
		syncQuESTEnv(env);
		if (env.rank==0) wtime_start = system_timer();	

		sol = solve3SAT(
			// inputs
			equ, numBools, numCycles, numClauses, thetaFactors,
			// outputs
			&totalSolProb, &runTime, &qubitZeroProbs, &clausePassProbs, &clauseReachedFailProbs,
			// progress monitoring
			progressFlag, progressFN);
	
		// end timing	
		syncQuESTEnv(env);
		if (env.rank==0) {
		    wtime_stop = system_timer ();
		    printf("\n\nTOTAL TIME = %.5f s\n\n", wtime_stop-wtime_start);
		    timingVec[trial]=wtime_stop-wtime_start;
		}
	}

	printf("local: %d, dist: %d\n", numLocal, numDist);
        // report timing to file
        if (env.rank==0){
                REAL totTime, avg, standardDev, temp, max, min;
                max=0; min=10e5;
		totTime=0;
		for (trial=0; trial<N_TRIALS; trial++){
			temp=timingVec[trial];
			totTime+=temp;
			if (temp<min) min=temp;
			if (temp>max) max=temp;
		}
		avg = totTime/(REAL)N_TRIALS;
		standardDev=0;
		for (trial=0; trial<N_TRIALS; trial++){
			temp = timingVec[trial]-avg;
			standardDev += temp*temp;
		}
		standardDev = sqrt(standardDev/(REAL)N_TRIALS);
		fprintf(timing, "%.8f, %.8f, %.8f, %.8f\n", avg, standardDev, max-avg, avg-min);
        }

	if (env.rank==0){
		// report state reached
		printf("reached ");
		for (int i=0; i < numBools; i++) 
			printf("%d", sol[i]);
		printf("> in %ds not including memory allocation ", runTime);
		printf("with probability %f\n", totalSolProb);
	
		if (isValidSolution(equ, sol, numClauses))
			printf("which was the correct state\n");
		else
			printf("WHICH WAS AN INCORRECT SOLUTION\n");
		
		// write input data to file
		createOrClearFile(outFN);
		saveInteger(
			numBools, 
			outFN, "numBools");
		saveInteger(
			numClauses, 
			outFN, "numClauses");
		saveInteger(
			numCycles, 
			outFN, "numCycles");
		saveDouble(
			initThetaFactor, 
			outFN, "initThetaFactor");
		saveInteger(
			updateTheta, 
			outFN, "updateThetaFlag");
		saveDoubleArray(
			thetaFactors, numCycles,
			outFN, "thetaFactors");
		
		// write output data to file
		saveInteger(
			runTime, 
			outFN, "runTime");
		saveDouble(
			totalSolProb, 
			outFN, "totalSolProb");
		saveIntegerArray(
			sol, numBools, 
			outFN, "quantumSol");
		saveDoubleArray(
			clauseReachedFailProbs, numClauses*numCycles,
			outFN, "clauseReachFailProbsFlat");
		saveDoubleArray(
			clausePassProbs, numClauses*numCycles,
			outFN, "clausePassProbsFlat");
		saveDoubleArray(
			qubitZeroProbs, numBools*numCycles,
			outFN, "qubitZeroProbsFlat");
	}	
	
	/*
	 * FREE MEMORY
	 */

	return EXIT_SUCCESS;

	free(equ);
	free(sol);
	free(inFN);
	free(outFN);
	free(progressFN);
	free(thetaFactors);
	free(qubitZeroProbs);
	free(clausePassProbs);
	free(clauseReachedFailProbs);
	
	closeQuESTEnv(env);

	return EXIT_SUCCESS;
}
