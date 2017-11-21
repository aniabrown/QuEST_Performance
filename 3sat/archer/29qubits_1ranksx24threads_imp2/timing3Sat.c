/** @file 
 * 3SAT simulation code
 */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <unistd.h>
# include <string.h>
# include <omp.h>

# include "../../../QuEST_v0.9.0/precision.h"
# include "../../../QuEST_v0.9.0/qubits.h"


# include <stdlib.h>
# include <sys/time.h>

//! Max number of qubits in the system
# define MAX_QUBITS   40
//! Max number of complete cycles that the user can reqest
# define MAX_CYCLES   300
//! Max number of complete clauses in the user-specified file defining the 3SAT
# define MAX_CLAUSES 300
//! 1: print minute details, 0: don't do that
# define STREAM_OF_CON 0

# define REPORT_TIMING 1
# define N_TRIALS 1

const long REAL Pi = 3.14159265358979323846264338327950288419716939937510;
const long REAL PiOver2 = 1.570796326794896619231321691639751442099;
const long REAL PiOver4 = 0.7853981633974483096156608458198757210493;
const long REAL oneOverRoot2= 0.7071067811865475244008443621048490392848;

// ==================================================================== //
//                                                                      //
//     system_timer -- precision walltime function, computes            //
//                     walltime based on the gettimeofday() function    //
//                                                                      //
// ==================================================================== //

REAL system_timer (void) {


        struct timeval time;

        gettimeofday (&time, NULL);

        return time.tv_sec + time.tv_usec / 1000000.0;

}


REAL interpolator(int mode, REAL initTheta, REAL progressPoint, int thisQ) {

	REAL theta=0;
	if (mode==0){
  		theta=initTheta;
  	}else{
  		if (mode==1){
	  		theta=initTheta+(PiOver2-initTheta)*progressPoint;
  		}else{
  			theta=initTheta+(PiOver2-initTheta)*progressPoint*progressPoint*progressPoint;
  		}
  	}
	return theta;
}

//--------------------------------------------------------------
//---------------------- START OF main()  ----------------------
//--------------------------------------------------------------
int main ( int argc, char *argv[] ){

	//
	// ===== INITIALISATION
	//
	
	// INIT ENVIRONMENT: ALWAYS REQUIRED ONCE AT BEGINNING OF PROGRAM
	// These two lines will automatically set up the environment (multinode,
	// openMP only etc)  
	QuESTEnv env;
	initQuESTEnv(&env);	
	
	//declare key varaibles

	REAL initTheta=0.5*PiOver2; //this is the angle theta for the algorithm, in radians of course, this default value is overwritten by user choice below
	int numCycles=50; //this sets how many complete cycles of operators we want to perform
	int orderedCycles=1; //this determines if we want to go through each cycle of clause checks in order (yes, when we set orderedCycles=1) the alternative is to choose each clause operator at random (orderedCycles=0) or random within cycles (orderedCycles=2). 
	REAL rError=0.00; //this is the amount of error to introduce on every rotation. NOT YET IMPLEMENTED
 
	REAL desiredFrac=0.3;
	
	int userSetRandom=0;
	int userSpecifiedFile=0;
//	int dumpAllProbs=0;
	int evolveThetaValue=0;
	int filterMode=0;
	long long int baseForPlots=-1;
	int fullMonitor=0;
	int printClauseProbs=0;
	
	int verboseWorking=STREAM_OF_CON;  //but user can switch verboseWorking to 1 so can be 1 even if S_O_C=0
	
	unsigned rSeed=time(0);
	unsigned long userRND=0;


	if (argc<2){
		printf("\nA usage example is './3SAT myFile -r=5 -s=100 -a=0.49 -verbose'\n");
		printf("where: myFile is the only mandatory argument and names the file containing the clause specification,\n");
		printf("  -verbose optionally specifies running commentary\n");
		printf("  -e is whether to Evolve theta value from initial (see -a) to pi/2 [default 0='no', 1 or 2='yes, linear or special']\n");
		printf("  -r is an adjustment to the random seed from clock [default 0]\n");
		printf("  -s is the number of full cycles of simulation [default 20]\n");
		printf("  -F is whether to use filter mode (one fewer qubit, fewer gates) 0='no', 1='yes' [default 1]\n");
		printf("  -a is the desired initial angle in the target state, as a portion of pi/2 [default 0.3]. \n");
		printf("  -p is whether to print the probs at per-clause level [default 0='no']. If p=N, N>1, the first N clause(s) are printed. \n");
		printf("  -d is whether to use the clauses in orDered sequence (d=1, the default), or to choose each clause randomly (d=0). \n");
		printf("  -M is whether to monitor the fidelity -v- solution, and the per-qubit basis, after each cycle [ 0(default)=no, 1=yes]. \n");
		exit(1);
	}
	
	int i; int skipArg=-1;
	
	printf("commandlineInvocation=\"");
	for (i=1; i<argc;  i++){
		printf(" %s",argv[i]);
	}
	printf("\";\n");
	
	for (i=1; i<argc;  i++)	if (strcmp(argv[i],"-verbose")==0)
	{
		printf("\n----------------------------------------\nSetting verbose mode active and starting run.\n");
		verboseWorking=1;
		skipArg=i;
	}
	for (i=1; i<argc;  i++){
		if ((argv[i][0]=='-') && (i!=skipArg)){
			int j=1;
			while (argv[i][j]=='-') j++;
			int failedToParse=1;
			char thisFlag=argv[i][j]; j++;
			if (argv[i][j]=='='){
				while (j>=0){
					argv[i][j]=' ';
					j--;
				}
				if (thisFlag=='a'){
					desiredFrac=strtof(argv[i], NULL);
					if ((desiredFrac<0) || (desiredFrac>=1)){	printf("\nSaw param a=%s. Target fraction of pi/2 must be between 0 and 1, exclusive. ABORTING.\n",argv[i]); exit(1);}
					if (verboseWorking!=0) printf("Command line argument specifies target fraction of pi/2 as %f.\n", desiredFrac);
					initTheta=desiredFrac*PiOver2;
					if (initTheta==0){ initTheta=0.00000001;}
					failedToParse=0;
				}
				if (thisFlag=='s'){
					int userCycles=atoi(argv[i]);
					if (userCycles>0){
						numCycles=userCycles;
						if (numCycles>MAX_CYCLES){ printf("\nTotal number of cycles requested %d exceeds the maxNumCycles. ABORTING.\n",numCycles); exit(1);}
						if (verboseWorking!=0) printf("Command line argument sets total cycles to %d.\n", userCycles);
						failedToParse=0;
					}else{
						printf("\nNumber of steps must be >= 1.ABORTING.\n");
						exit(1);
					}
				}
				if (thisFlag=='p'){
					printClauseProbs=atoi(argv[i]);
					if (printClauseProbs>=0){
						if (verboseWorking!=0) printf("Command line argument set p=%d relating to per-clause recording.\n", printClauseProbs);
						failedToParse=0;
					}else{
						printf("\nSetting for -p must be '0' (mode off) or '1' (mode on for all), or >1 (mode on for N clauses).ABORTING.\n");
						exit(1);
					}
				}
				if (thisFlag=='F'){
					filterMode=atoi(argv[i]);
					if ((filterMode==0) || (filterMode==1)){
						if (verboseWorking!=0) printf("Command line argument filter mode status to %d.\n", filterMode);
						failedToParse=0;
					}else{
						printf("\nSetting for filter mode must be '0' (mode off) or '1' (mode on).ABORTING.\n");
						exit(1);
					}
				}
				if (thisFlag=='e'){
					evolveThetaValue=atoi(argv[i]);
					if ((evolveThetaValue>=0) && (evolveThetaValue<=2)){
						if (verboseWorking!=0) printf("Command line argument sets evolution switch to %d.\n", evolveThetaValue);
						failedToParse=0;
					}else{
						printf("\nIf -e=x is specified, x must be 0,1 or 2. ABORTING.\n");
						exit(1);
					}
				}
				if (thisFlag=='M'){
					fullMonitor=atoi(argv[i]);
					if ((fullMonitor>=0) && (fullMonitor<=1)){
						if (verboseWorking!=0){ 
							printf("Command line argument specifies whether to monitor: "); 
							fullMonitor?printf("yes\n"):printf("no\n");
						}
						failedToParse=0;
					}else{
						printf("\nIf -M=x is specified, x must be 0 or 1. ABORTING.\n");
						exit(1);
					}
				}
				if (thisFlag=='r'){
					userRND=strtoul(argv[i], NULL, 0);
					if (verboseWorking!=0) printf("Command line argument adjusts the clock-based RND seed by +100*%lu.\n", userRND);
					rSeed+=100*userRND;
					userSetRandom=1;
					failedToParse=0;
				}
				if (thisFlag=='n'){
					REAL noiseLevel=strtof(argv[i], NULL);
					if ((noiseLevel<0) || (noiseLevel>10)){	printf("\nSaw param n=%s. Noise level must be in range 0 to 10 percent inclusive. ABORTING.\n",argv[i]); exit(1);}
					if (verboseWorking!=0) printf("Command line argument specifies noise level as %f percent.\n", noiseLevel);
					rError=noiseLevel/100.; //percent to absolute proportion
					failedToParse=0;
				}
				if (thisFlag=='d'){
					orderedCycles=atoi(argv[i]);
					if ((orderedCycles<0) || (orderedCycles>2)){	printf("\nIf -d=x is specified, then x must be 0, 1 or 2.\n"); exit(1);}
					if (verboseWorking!=0) printf("Command line argument '-d=%d' has been specified to determine where clause sequence is ordered or random.\n",orderedCycles);
					failedToParse=0;
				}
//				if (thisFlag=='m'){
//					cyclesToTrackAlign=atoi(argv[i]);
//					if (cyclesToTrackAlign<0){	printf("\nIf -d=x is specified, then x must be 0 [default, meaning 'off'] or positive.\n"); exit(1);}
//					if (verboseWorking!=0) printf("Command line argument '-m=%d' has been specified to how many cycles to track at the per-qubit level.\n",cyclesToTrackAlign);
//					failedToParse=0;
//				}
				if (thisFlag=='b'){
					baseForPlots=atoll(argv[i]);
					if (baseForPlots<0){	printf("\nIf -b=n is specified, then x must be 0 or positive; by default b is not used.\n"); exit(1);}
					if (verboseWorking!=0) printf("Command line argument '-b=%Ld' has been specified, forcing per-qubit plots to be w.r.t. soln %Ld.\n",baseForPlots,baseForPlots);
					failedToParse=0;
				}
			}
			if (failedToParse==1){ printf("\nCannot interpret command line argument %s. (Note use '=' as in '-s=100'). ABORTING.\n",argv[i]); exit(1); }
		
		}else{ // an argument not starting with - detected
			if (i!=skipArg){
				if (userSpecifiedFile!=0){ printf("\nTwo file names specified ? I see %s and %s. ABORTING.\n",argv[userSpecifiedFile],argv[i]); exit(1); }
				userSpecifiedFile=i;
			}	
		}
	}
	
	
    srand(rSeed);
	if (verboseWorking!=0) printf("Random seed=%i\n",rSeed);
	
	
	if (verboseWorking!=0) printf("Desired target angle=%f(pi/2)\n",desiredFrac);
	
  
	
	int numClauses=0;
	int numBits=0;
	
	int indexList[MAX_CLAUSES][3]={{0}};      //indexList is an array that will contain the same information as "clauses" but without the minus signs
	int notOrDirect[MAX_CLAUSES][3]={{0}};    //notOrDirect is an array that will have an entry for each clause and each bit within the clause, it will be 0 if the clause requires the bit to be TRUE and 1 if it wants it to be FALSE. So it is storing the sign information from the array 'clauses'

	
	
//---- the following block of code reads in the clauses specifying the 3SAT from the user's named file	

	if (userSpecifiedFile==0){
		printf("Command line does not specify a file containing the clauses; aborting.\n"); exit(1);
	}else{
		char thisFileName[200];
		strcpy(thisFileName,argv[userSpecifiedFile]); 
		if (verboseWorking!=0) printf("The file containing clauses is named %s.\n", thisFileName); 
		FILE * fPtr=fopen(thisFileName,"r");
		if (fPtr==NULL){ printf("Error openning the file named %s.\n", thisFileName); exit(1); }
		char newLine[30];
		char oldLine[30];
		int clCount=0;
		while(!feof(fPtr)){	
			fgets(newLine,30,fPtr);
			int len=strlen(newLine);
			if ((len>1) && (strcmp(newLine,oldLine) != 0)){
				char * pEnd;
				int a = (int) strtol ( newLine, &pEnd,10);
				int b = (int) strtol ( pEnd, &pEnd,10);
				int c = (int) strtol ( pEnd, &pEnd,10);
				
				if (STREAM_OF_CON) printf("Read in clause %d %d %d\n",a,b,c);
				
				if ((a==0) || (b==0) || (c==0)){
					printf("Could not read three non-zero values when attempting to read clause %d from line \"%s\". ABORTING.\n",clCount,newLine); exit(1);
				}
				// now convert clauses into distinct arrays for index (c style, 0..n-1) and NOT status
				// also make the order [0],[1],[2] be greatest index to smallest
				
				notOrDirect[clCount][0] = 0;
				notOrDirect[clCount][1] = 0;
				notOrDirect[clCount][2] = 0;
				if (a<0){notOrDirect[clCount][0] = 1; a=-a;}
				if (b<0){notOrDirect[clCount][1] = 1; b=-b;}
				if (c<0){notOrDirect[clCount][2] = 1; c=-c;}
				
				
				if ((a==b) || (b==c) || (a==c) ){
					printf("Repeated variable in clause %s. ABORTING.\n",newLine); exit(1);
				}
				
				if (a>numBits) numBits=a;
				if (b>numBits) numBits=b;
				if (c>numBits) numBits=c;
				
				indexList[clCount][0] =a-1; 
				indexList[clCount][1] =b-1; 
				indexList[clCount][2] =c-1;
				
				clCount++;
				
				if (clCount>MAX_CLAUSES){ printf("Number of clauses %d exceeds the maximum. ABORTING.\n",clCount-1); exit(1);}
			}
			strcpy(oldLine,newLine);
		}
		numClauses=clCount;
		if (verboseWorking!=0) printf("Reached end of file; read %d clauses.\n",numClauses);
		fclose(fPtr);
	}

//end of clause file reader -------------------------------------

	
	

//--------- the following code block seeks the classical solution to the 3SAT. Needed only so as to plot the fidelity metrics
	
	if (verboseWorking) printf("Now seeking solution(s) by classical brute force.\n");

	long long int binCases = 1LL << numBits; //this simply computes 2^numBits 
	int numSolutions=0;
	long long int oneSolution;
	long long int index;
	# ifdef _OPENMP
		# pragma omp parallel for \
		default  (none) \
		shared   (indexList,notOrDirect,binCases,numClauses,numSolutions,oneSolution,numBits,verboseWorking	) \
		private  (index) \
		schedule (static)
	# endif 
	for (index=0; index<binCases; index++){ 
		int bitString[MAX_QUBITS]; //here we define an array which will store the binary bits in the number 'index'
  		//note, using fixed array def above seems to yield much faster performance than malloc(numQbits, ..)
  		long long int tmp= index;
  		int sweepBits;
		for (sweepBits=0; sweepBits<numBits; sweepBits++){
			bitString[sweepBits]= tmp & 1LL;
			tmp = tmp >> 1;
		}
		int satisfies=1;
		int thisCl;
		for (thisCl=0; (thisCl<numClauses) && (satisfies==1); thisCl++){ //this nested loop will detect if the current bit string is a solution to the 3SAT
			int valid=0;
			int thisEle;
			for (thisEle=0; (thisEle<3) && (valid==0); thisEle++){
				int thisBit=indexList[thisCl][thisEle];
				int feared=notOrDirect[thisCl][thisEle];
				if (bitString[thisBit]!=feared){ 
					valid=1; //nb we use [thisBit-1] since the 1st bit is slot [0]
				} 
			}
			if (valid==0){
				satisfies=0; 
			}
		}
//  	printf("trying clause %d and finding satisfies=%d\n",thisCl,satisfies);
		if (satisfies==1){
  			oneSolution=index;
  			numSolutions++;
  			if (verboseWorking==1) printf("FOUND solution i.e. %Ld\n",oneSolution);
  		 	if (STREAM_OF_CON){
  		 		printf("bin: ");
  		 		for (sweepBits=0; sweepBits<numBits; sweepBits++) printf("%d",bitString[sweepBits]);
  		 		printf("\n");
  		 	}
  		}
	}//end soln search loop
	
	
	if (numSolutions>1){ printf("Multiple Solutions! %d were found; one example is %Ld.\n",numSolutions,oneSolution); exit(1);} //if there is more than one solution, or zero solutions, then this will abort the program (because what happens next assumes exactly one)
	if (numSolutions==0){ 
		if (baseForPlots>-1){
			printf("alert=\"NO CLASSICAL SOLUTION EXISTS! Proceeding with user-specified base solution %Ld.\";\n",baseForPlots);
			oneSolution=baseForPlots;
		}else{
			printf("alert=\"NO CLASSICAL SOLUTION EXISTS! Proceeding as if soln=000..0.\";\n");
			oneSolution=0;
		}
	}else{
		if (baseForPlots>-1){
			printf("alert=\"WARNING! Overriding the true solution %Ld with the user-specified base solution %Ld.\";\n",oneSolution,baseForPlots);
			oneSolution=baseForPlots;
		}
	}
	
	int * solnString = (int*) calloc(numBits, sizeof(int)); //here we define an array which will store the binary bits in the number 'index'
  	long long int tmp= oneSolution;
  	int sweepBits;
	for (sweepBits=0; sweepBits<numBits; sweepBits++){
		solnString[sweepBits]= tmp & 1LL;
		tmp = tmp >> 1;
	}
  
	if  (verboseWorking!=0){
		printf("Sole classical 3SAT solution is %Ld which is ",oneSolution);
		int cntDwn;
		for (cntDwn=numBits-1; cntDwn>=0; cntDwn--) printf("%d",solnString[cntDwn]); //here the counter goes down from numQubits-2 because it is easier for the user to understand the binary number if he sees printed it with the 0th bit on the left, then the 1st, the 2nd, etc
	  	printf("\n");
  	}
  	
  	//so now we have found the classical solution and stored it in the array called 'solnString'
	    REAL wtime_start, wtime_stop;
/*
	    REAL *timingVec;
	    FILE *timing;
	    int trial;
	    char envString[255];
	    char filename[255];
	    if (REPORT_TIMING && env.rank==0) timingVec = (REAL*) malloc(N_TRIALS*sizeof(timingVec));

	    if (REPORT_TIMING && env.rank==0){
		MultiQubit testMQ;
		createMultiQubit(&testMQ, numBits, env);
		getEnvironmentString(env, testMQ, envString);
		destroyMultiQubit(testMQ, env);
		sprintf(filename, "TIMING3SAT_%s.csv", envString);
		timing = fopen(filename, "w");
		fprintf(timing, "qubit, time(s), standardDev, maxDelta, minDelta\n");
	    }
*/

	// --------------------- BEGIN TIMING
	if (REPORT_TIMING) syncQuESTEnv(env);
	if (REPORT_TIMING && env.rank==0) wtime_start = system_timer();
	
	int numQubits=numBits;
//	long long int numAmps = 1LL << numQubits; //this simply computes 2^numQubits and puts that number into the variable numAmps (short for number of amplitudes)



//	REAL * fidStore = (REAL *) calloc(numCycles, sizeof(REAL));
//	REAL * TStore = (REAL *) calloc(numCycles, sizeof(REAL));


	REAL * bitRec = (REAL *) calloc(numCycles*numBits, sizeof(REAL));
	
	REAL * perClauseFail= (REAL *) calloc(numCycles*numClauses, sizeof(REAL));
	
	REAL * currentOrientation= (REAL *) calloc(numQubits, sizeof(REAL));
	REAL * targThetaPerCycPerQ= (REAL *) calloc(numCycles*numQubits, sizeof(REAL));

	

	// CREATE QUBIT OBJECT: REQUIRED ONCE PER MULTIQUBIT OBJECT	
	// Before doing any operations on a set of qubits, create the MultiQubit object that will be used to 
	// represent the qubits
	MultiQubit qubitRegister; 
	printf("Number of qubits: %d\n", numQubits);
	createMultiQubit(&qubitRegister, numQubits, env);
	
	// Reporting
	if  (verboseWorking!=0 && env.rank==0){
		printf("-------------------\nInstruction \"reportMultiQubitParams(qubitRegister)\" yields:\n");
		reportMultiQubitParams(qubitRegister);
		printf("Instruction \"reportQuESTEnv(qubitRegister)\" yields:\n");
		reportQuESTEnv(env);
		printf("-------------------\n");
	}
	
	// initialise the state to |0000..0>
	initStateZero(&qubitRegister);   
	
	if (env.rank==0) if (verboseWorking!=0) printf("Completed initStateVec\n");
	
	Complex alpha, beta;
  	alpha.real=0; alpha.imag=oneOverRoot2; beta.real=0; beta.imag=oneOverRoot2; //H rotation, Hadamard-like
	int thisQ;
	for (thisQ=0; thisQ<numQubits; thisQ++){  //this is used if we want initial state |++..+>
//		rotateQubit(qubitRegister, thisQ, alpha, beta);
//		if (verboseWorking!=0) printf("Completed rotation to |+> for Q%d\n",thisQ);
	}
	
//	REAL lowestBR=1.; 
//	int lowPointBR=-1;
//	int noneUnder51=0;
//  int cycForAllOver51=-1;
	  
	int thisCycle;
	for (thisCycle=0; thisCycle<numCycles; thisCycle++){  	///loop to populate targTheta
		for (thisQ=0; thisQ<numQubits; thisQ++){  			
			targThetaPerCycPerQ[thisCycle*numQubits+thisQ]=interpolator(evolveThetaValue,initTheta, ((REAL)thisCycle)/(numCycles-1),thisQ);
			if (STREAM_OF_CON) printf("%d:%d:(%d):%f ",thisCycle,thisQ,thisCycle*numQubits+thisQ,targThetaPerCycPerQ[thisCycle*numQubits+thisQ]);
		}
	}
	
  		// ********** *********** *********** **********
  		// ********** starting Q algorithm ***********
  		// ********** *********** *********** **********
	
 

  	long int totClausesDone=0;
  	REAL oneRunSuccProb=1.;
	for (thisCycle=0; thisCycle<numCycles; thisCycle++){  ///start of CYCLE
 		if  (verboseWorking!=0 && env.rank==0){
 			printf("----------------- Starting cycle %d --------------------\n",thisCycle);
			printf("curOrient, targTheta for each qubit are: \n");
			for (thisQ=0; thisQ<numQubits; thisQ++) printf("Q%d (%d) : %f, %f; ",thisQ,thisCycle*numQubits+thisQ,currentOrientation[thisQ],targThetaPerCycPerQ[thisCycle*numQubits+thisQ]);
   			printf("\n");
		}
  		
  		// ********** *********** *********** **********
  		// ********** starting a clause cycle  ***********
  		// ********** *********** *********** **********
  		
	  	REAL cycleSuccProb=1.;
	  	int thisCl;
  		for (thisCl=0; thisCl<numClauses; thisCl++){
			
			if (env.rank==0) if (verboseWorking==1) printf(">> Performing clause %d.\n",thisCl);
			
			int thisEle; 
			
			for (thisEle=0; thisEle<3; thisEle++){ //in this loop, we go through each of the three element (bits) in the clause and rotate the qubits accordingly. 
				int thisQ=indexList[thisCl][thisEle];
				if (env.rank==0) if (STREAM_OF_CON) printf("Using element %d whose notOrDir state is %d\n", thisQ, notOrDirect[thisCl][thisEle]);
				REAL effectiveTheta=0;
				REAL theta=(1-2*notOrDirect[thisCl][thisEle]); //just sets +1 / -1
				theta = theta*targThetaPerCycPerQ[thisCycle*numQubits+thisQ];
				
//				theta += PiOver2; //needed if initial state |+++..+> as we can't efficiently filter out |---> so we will filter |111> presently; we need not rotate

				effectiveTheta = +theta-currentOrientation[thisQ];
				
				if ((effectiveTheta!=0.) && (fabs(effectiveTheta) < 0.00000001)){ 
					printf("NEGLIGABLE THETA=%f POSSIBLE ERROR\n",effectiveTheta); 
					exit(1);
				}
				currentOrientation[thisQ] = theta;
				
				if (effectiveTheta!=0.){
					if (env.rank==0) if (verboseWorking!=0) printf("Rotating Q%d;\n",thisQ);
  		
  					Complex alpha, beta;
  	
					alpha.real=cos(effectiveTheta/2); alpha.imag=0; beta.real=-sin(effectiveTheta/2); beta.imag=0; //specifying Y rotation
				
					rotateQubit(qubitRegister, thisQ, alpha, beta);
				} else if (verboseWorking!=0 && env.rank==0) printf("No rotation needed for Q%d;\n",thisQ);
			}//end of rotation of three qubits
			
			if (verboseWorking!=0 && env.rank==0) printf("Now filter111 on Q%d, Q%d, Q%d.\n",indexList[thisCl][0], indexList[thisCl][1] ,indexList[thisCl][2]);
			
			REAL probOfZero=filterOut111(qubitRegister, indexList[thisCl][0], indexList[thisCl][1] ,indexList[thisCl][2]);


			cycleSuccProb=cycleSuccProb*probOfZero;
			perClauseFail[totClausesDone]=1-probOfZero; 
			totClausesDone++;
			
		   if (STREAM_OF_CON && env.rank==0) printf("Clause %d probOfZero=%f thus 1-p=%f\n",thisCl,probOfZero,1-probOfZero);
		   
		   if (0){      //per CLAUSE qubit orientation tracking!
			int sweepQ;
			if (verboseWorking!=0 && env.rank==0) printf("Now rotating to basis for alignment check:\n");
			for (sweepQ=0; sweepQ<numQubits; sweepQ++){
				
//				REAL effectiveTheta = -PiOver4 - currentOrientation[sweepQ]; //this for initial |+++..+>
				
				REAL effectiveTheta = -PiOver2 - currentOrientation[sweepQ];
				
				if ((effectiveTheta!=0.) && (fabs(effectiveTheta) < 0.00000001)){ 
					printf("NEGLIGABLE THETA=%f POSSIBLE ERROR\n",effectiveTheta); 
					exit(1);
				}
				currentOrientation[sweepQ] = -PiOver2;
				
				if (effectiveTheta!=0.){
					if (STREAM_OF_CON && env.rank==0) printf("Calling for rotation of Q%d; ",sweepQ);
  		
  					Complex alpha, beta;
  	
					alpha.real=cos(effectiveTheta/2); alpha.imag=0; beta.real=-sin(effectiveTheta/2); beta.imag=0; //specifying Y rotation
				
					rotateQubit(qubitRegister, sweepQ, alpha, beta);
				} else if (STREAM_OF_CON && env.rank==0) printf("No rotation needed for qubit %d\n",sweepQ);
				
				if (STREAM_OF_CON && env.rank==0) printf("Finding prob of zero... for Q%d\n",sweepQ);
				printf("%f",findProbabilityOfOutcome(qubitRegister, sweepQ, 0));
				if (sweepQ<numQubits-1){ printf(","); }else{ if (thisCycle<numCycles-1) printf("},{"); else printf("}};");}
			}
			}
			if (verboseWorking!=0 && env.rank==0) printf("-----------------\n");
		} //ending clause cycle

		if (fullMonitor){
			int sweepQ;
			if (verboseWorking!=0 && env.rank==0) printf("Now rotating to basis for alignment check:\n");
			for (sweepQ=0; sweepQ<numQubits; sweepQ++){
				
//				REAL effectiveTheta = -PiOver4 - currentOrientation[sweepQ]; //this for initial |+++..+>
				
				REAL effectiveTheta = -PiOver2 - currentOrientation[sweepQ];
				
				if ((effectiveTheta!=0.) && (fabs(effectiveTheta) < 0.00000001)){ 
					printf("NEGLIGABLE THETA=%f POSSIBLE ERROR\n",effectiveTheta); 
					exit(1);
				}
				currentOrientation[sweepQ] = -PiOver2;
				
				if (effectiveTheta!=0.){
					if (STREAM_OF_CON && env.rank==0) printf("Calling for rotation of qubit %d\n",sweepQ);
  		
  					Complex alpha, beta;
  	
					alpha.real=cos(effectiveTheta/2); alpha.imag=0; beta.real=-sin(effectiveTheta/2); beta.imag=0; //specifying Y rotation
				
					rotateQubit(qubitRegister, sweepQ, alpha, beta);
				} else if (STREAM_OF_CON && env.rank==0) printf("No rotation needed for qubit %d\n",sweepQ);
				
				if (verboseWorking!=0 && env.rank==0) printf("Finding prob of zero for Q%d to be: ",sweepQ);
				bitRec[thisCycle*numQubits+sweepQ] = findProbabilityOfOutcome(qubitRegister, sweepQ, 0);
				if (verboseWorking!=0 && env.rank==0) printf("%f\n",bitRec[thisCycle*numQubits+sweepQ]);
			}
  			//rotate and compute fidelity
  		}
  		
		oneRunSuccProb=oneRunSuccProb*cycleSuccProb;
  	
	}
	if (REPORT_TIMING) syncQuESTEnv(env);	
        if (REPORT_TIMING && env.rank==0) { 
            wtime_stop = system_timer ();
	    printf("\n\nTOTAL TIME = %.5f s\n\n", wtime_stop-wtime_start); 
        }


	if (env.rank==0) printf("oneRunSuccProb=%f\n",oneRunSuccProb);

	if (fullMonitor!=0 && env.rank==0){
		printf("align[%ld]={{",userRND);
		int toTrack;
		for (toTrack=0; toTrack<numCycles; toTrack++) {
			for (sweepBits=0; sweepBits<numQubits; sweepBits++) {
				printf("%f", bitRec[toTrack*numQubits+sweepBits]);
				if (sweepBits<numQubits-1) printf(",");
	  		}
  			if (toTrack<numCycles-1) printf("},{");
		}
		printf("}};\n");
	}

	if (env.rank==0){	
		if (printClauseProbs) printf("failProbOverClauses[%ld]={",userRND);
		int thisCl;
		REAL rollingSucc=1.;
		REAL * cumulat = (REAL *) calloc(numCycles, sizeof(REAL));
		REAL avClCount=0;
		long int count=0;
		for (thisCl=0; thisCl<numCycles*numClauses; thisCl++) {
				
			REAL thisVal = rollingSucc*perClauseFail[thisCl];
			avClCount += (thisCl+1)*thisVal;	//+1 because if a clause fails, we still did the work for that clause
			rollingSucc=rollingSucc*(1-perClauseFail[thisCl]);
			if ((printClauseProbs==1) || (printClauseProbs>count) ){
				if ((thisVal>0.0001) || (thisVal<=0)){ //clearly should not be <0 tho!
					printf("%f", thisVal);
				}else{
					printf("%.4e", thisVal);
				}
				if (thisCl<numCycles*numClauses-1) printf(",");
			}
			count++;
			if (count % numClauses ==0)
				cumulat[(count/numClauses)-1]=rollingSucc;
			
		}
		if (printClauseProbs) printf("};\n");
		cumulat[numCycles-1]=rollingSucc;
		
		printf("cumulatSuccPerCyc[%ld]={",userRND);
		for (thisCycle=0; thisCycle<numCycles; thisCycle++) {
			REAL thisVal=cumulat[thisCycle];
			if ((thisVal>0.001) || (thisVal<=0)){ //clearly should not be <0 tho!
				printf("%f", thisVal);
			}else{
				printf("%.4e", thisVal);
			}
			if (thisCycle<numCycles-1) printf(",");
		}
		printf("};\n");
		
		printf("runSuccessProb[%ld]=%f\n",userRND,rollingSucc);
		avClCount += numCycles*numClauses*rollingSucc; //accounts for the portion of cases that succeed to end
		printf("avClCount[%ld]=%f\n",userRND,avClCount/rollingSucc);
	}
	
	
	//
	// ======== CLEANUP
	//
	
	// free memory

	// REQUIRED ONCE PER MULTIQUBIT OBJECT	
	// When all operations on a set of qubits are completed, destroy the object
	destroyMultiQubit(qubitRegister, env);


	// ALWAYS REQUIRED ONCE AT END OF PROGRAM: 
	// These two lines will perform any necessary cleanup of the environment (multinode,
	// openMP only etc)  
	closeQuESTEnv(env);

	return EXIT_SUCCESS;

}

