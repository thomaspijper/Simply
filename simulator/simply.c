/* 
 *  Copyright (c) 2005 Gabriele Keller
 *  Copyright (c) 2006 Don Stewart
 *  Copyright (c) 2007 Hugh Chaffey-Millar
 *  Copyright (c) 2016 Thomas Pijper
 *
 *  This file is part of Simply.
 *
 *  Simply is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  Lesser General Public License for more details.
 */

// Headers usable under both Windows and Linux
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <inttypes.h>
#include "genpolymer.h"
#include "simply.h"
#include <assert.h>
#include <limits.h>
#include <float.h>
#include "dSFMT.h"

// Debugging functions
//#define DEBUG 1
//#define TRACE 1
//#define NO_COMM 1
#define MONO_AUDIT 1
#define EXPLICIT_SYSTEM_STATE 1

// Control of random number generation
//#define CHANGE_SEED 1

// Control of system time calculation
#define MICRO_CALC_TIME 1

//#define SCALING 1

#define VERBOSE 0
#define IFVERBOSE if (VERBOSE == 1)
#define RANK if (myid == 0)

// Aliases used for printing to files
#define FULL 0
#define PROFILES 1
#define START 2

#define HOST_ID 0
#define START_MWD_SIZE 512 // must be a power of 2
#define INIT_STATE_COMM_SIZE (6 * sizeof(pcount) * START_MWD_SIZE)
#define TREE_SAMPLE	10
#define PACKET_UPSCALE_FACTOR 1.7F // factor by which to increase packet size if it is too small
#define PACKET_SIZE_FACTOR 1.7F
#define MAX_FILENAME_LEN 100
#define MAX_FILE_SIZE 1048576

#define AVOGADRO           6.022140857E23

#define NUM_RANDS 10000 // amount of pseudorandom numbers to generate at a time; should not be smaller than 382 and must be an even number

// Messages for setting up data at start
#define SETUP_END 0
#define SETUP_CONVDATA 1
#define SETUP_SYSTEMSCALES 2

enum bools {False=0,True=1};

int numprocs = -1;
int64_t lastReduceTime = 0;
int currentComparisonComplexity = -1;

int checkPointing = 0;

static sysState state;


#if defined(__GNUC__)
// Defines some includes for GCC
#include <mpi.h>
#include <unistd.h>
#include <sys/time.h>
#include <getopt.h>
#include <sched.h>
//#include <mcheck.h>
#define _GNU_SOURCE
#define _BSD_SOURCE // required by timing functions
#elif defined(_MSC_VER)
#include <process.h>
// Specify alternatives to POSIX unique headers/functions
#include "mpi.h"
// Define gettimeofday() for Windows systems
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
int gettimeofday(timeval *tp, int *null_value) {
	static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

	SYSTEMTIME  system_time;
	FILETIME    file_time;
	uint64_t    time;

	GetSystemTime(&system_time);
	SystemTimeToFileTime(&system_time, &file_time);
	time = ((uint64_t)file_time.dwLowDateTime);
	time += ((uint64_t)file_time.dwHighDateTime) << 32;

	tp->tv_sec = (int64_t)((time - EPOCH) / 10000000L);
	tp->tv_usec = (int64_t)(system_time.wMilliseconds * 1000);
	return 0;
}
#if _MSC_VER < 1900
#error Visual Studio versions older than VS 2015 are not supported.
#endif
#else
#error The compiler you are trying to use is not supported.
#endif

// Define inline directives for MSVC and GCC
#if defined(_MSC_VER)
#define INLINE __inline
#elif defined(__GNUC__)
#define INLINE inline
#endif

// Declare fetchpid() for Windows and Linux systems
int fetchpid(void) {
#if defined(_MSC_VER)
	return _getpid();
#elif defined(__GNUC__)
	return getpid();
#endif
}

#if defined(__GNUC__)
// ****** the following stuff required for getopt on UNIX systems
int getopt(int argc, char * const argv[],const char *optstring);

pcount init_monomer_particles = NO_OF_MONOMER_MOLS;

extern char *optarg;
extern int optind, opterr, optopt;

int getopt_long(int argc, char * const argv[],
           const char *optstring,
           const struct option *longopts, int *longindex);

int getopt_long_only(int argc, char * const argv[],
           const char *optstring,
           const struct option *longopts, int *longindex);

// **************************** end of getopt *****
#endif

int64_t total_wtime = 0, total_rtime = 0, reduces = 0;

// Stores the header of the blocks of data that are communicated.
typedef struct {
    int 			stateTooBig;		// Flag indicating that size of 
										//  state has become too large.
    ptime           time;				// Time of system
	double			deltatemp;			// Temperature deviation of the system
	int 			noMoreReactions;	// Counts number of nodes that 
										//  have no events possible.
	pcount 			globalAllMonomer;	// Counts number of monomer 
										//  particles that have been 
										//  consumed and that still exist 
										//  as monomer.
} StatePacket;

const int startsize = START_MWD_SIZE;

int myid = -1;

void print_kinetic_model(void);
int monomerAudit(const char *str);

void print_state(void);
pcount check_tree(int no_of_leaves, int levelSize, int offset, pcount *mwd, int spec_ind);
       
INLINE int max_value(int x, int y) {
	return (x < y ? y : x);
}

INLINE pcount min_value(pcount x, pcount y) {
	return (x > y ? y : x);
}

// Variables needed for pseudorandom number generation
static w128_t dummy0[NUM_RANDS / 2 + 1];
static w128_t dummy1[NUM_RANDS / 2 + 1];
static w128_t dummy2[NUM_RANDS / 2 + 1];
static double *rndArray0 = (double *)dummy0;
static double *rndArray1 = (double *)dummy1;
static double *rndArray2 = (double *)dummy2;
static int rndCounter0 = NUM_RANDS;
static int rndCounter1 = NUM_RANDS;
static int rndCounter2 = NUM_RANDS;

/* Returns a pseudorandom number.
 *   randomProb(0) returns number in interval (0,1)
 *   randomProb(1) returns number in interval [0,1)
 *   randomProb(2) returns number in interval (0,1]
 *   randomProb(3) returns number in interval [0,1] -- currently (0,1)
 */
INLINE static double randomProb(int x) {

	if ((x == 0)
		|| (x == 3)) { // interval [0,1] is not available with dSFMT.h
		if (rndCounter0 == NUM_RANDS) {
			dsfmt_gv_fill_array_open_open(rndArray0, NUM_RANDS); // interval (0,1)
			rndCounter0 = 0;
		}
		double rnd = rndArray0[rndCounter0];
		rndCounter0 += 1;
		return rnd;
	}
	else if (x == 1) {
		if (rndCounter1 == NUM_RANDS) {
			dsfmt_gv_fill_array_close_open(rndArray1, NUM_RANDS); // interval [0,1)
			rndCounter1 = 0;
		}
		double rnd = rndArray1[rndCounter1];
		rndCounter1 += 1;
		return rnd;
	}
	else if (x == 2) {
		if (rndCounter2 == NUM_RANDS) {
			dsfmt_gv_fill_array_open_close(rndArray2, NUM_RANDS); // interval (0,1]
			rndCounter2 = 0;
		}
		double rnd = rndArray2[rndCounter2];
		rndCounter2 += 1;
		return rnd;
	}
	else {
		printf("You've reached an unavailable option for randomProb. Exiting...\n");
		exit(EXIT_FAILURE);
	}
}

/* Makes each process take the correct number of particles, ensuring conservation of
 * particles accross all systems.
 */
INLINE pcount takeSome(pcount total) {
		pcount ans = total/numprocs;
		pcount leftovers = total % (pcount)numprocs;
		if (myid < leftovers)
			ans++;
		return ans;
}

INLINE pcount leaveSome(pcount total) {
		pcount rough = total/numprocs;
		pcount leftovers = total % (pcount)numprocs;
		pcount correction = (myid >= leftovers) ? leftovers : myid;
		pcount ans = myid*rough + correction;
		return ans;
}

/* Deals with overflow of array containing mwd tree
 */
void double_arraysize (pcount **arr, int curr_size) {
  int newsize = 2 * curr_size;
  pcount *new_arr = calloc (sizeof (pcount) * (2 * newsize - 1), 1);
  pcount *old_arr = *arr;

  int p      = 1;  /* size of current level */
  int src_ix = 0;  /* src index             */
  int tgt_ix = 1;  /* target index          */

  new_arr[0] = old_arr[0];
  while (p < newsize) {
    memcpy (&(new_arr[tgt_ix]), &(old_arr[src_ix]), p * sizeof (pcount));
    src_ix += p;
    p = 2*p;
    tgt_ix = 2 * p -1;
  }
  *arr = new_arr;
  free (old_arr);
}

void startTimer(Timer *t) {
	gettimeofday(&t->start_t,NULL);
}

uint64_t readTimer(Timer *t) {
	gettimeofday(&t->end_t,NULL);
	return (1000000*t->end_t.tv_sec+t->end_t.tv_usec)-(1000000*t->start_t.tv_sec+t->start_t.tv_usec);
}

double readTimerSec(Timer *t) {
	gettimeofday(&t->end_t,NULL);
	double begin = (double)t->start_t.tv_sec+t->start_t.tv_usec/1000000.0;
	double end = (double)t->end_t.tv_sec+t->end_t.tv_usec/1000000.0;
	return (end-begin);
}

void dumpTree(mwdStore *x) {
	for (int i=1; i <= x->maxEntries; i*=2) {
		for (int j=0; j < i; j++) {
			printf (" %lld",x->mwd_tree[i-1+j]);
		}
		printf("\n");
	}
}

void dumpAllTrees() {
	for (int i =0; i < NO_OF_MOLSPECS; i++) {
		printf("%s:\n",name(i));
		for (int a=0; a < state.arms[i]; a++) {
			dumpTree(&state.mwds[i][a]);
		}
	}
}

INLINE void adjustTree(int spec_ind, mwdStore  *mwd, int leave_ind, int diff) {

	int ind;
	int levelSize = 1;

	while (leave_ind >= mwd->maxEntries) {
		printf("tree data structure overflow: %s leaveIx = %d maxEntries = %d\n", name(spec_ind),leave_ind, mwd->maxEntries);
		double_arraysize(&(mwd->mwd_tree), mwd->maxEntries);
		mwd->maxEntries *= 2;
	}

	while (levelSize <= mwd->maxEntries) {
		ind = (levelSize - 1) + (levelSize*leave_ind)/mwd->maxEntries;
		mwd->mwd_tree[ind] += diff;
		levelSize *= 2;
	}

#if defined(DEBUG)
	//	printf("adjustMolCnt(%d, %d, %d)...\n",spec_ind,leave_ind,diff);
	//	printf("check_tree running on %s\n", name(spec_ind));
	//	dumpTree(&state.mwds[spec_ind]);
//	check_tree(mwd->maxEntries, 1, 0, mwd->mwd_tree, spec_ind);
#endif
}

/* Add `diff' molecules of chain length (leave_ind + 1) of type
 * spec_ind into tree. Increase tree size if necessary. 
 *
 * From CodeGen polysim.c.
 */
INLINE void adjustMolCnt(int spec_ind, chainLen *lengths, int arms, int diff) {
	if (spec_ind < MAXPOLY) {
		adjustTree(spec_ind, &state.mwds[spec_ind][0], lengths[0]-1,diff);
	}
	else {
		for (int a=0; a<arms; a++) {
			adjustTree(spec_ind, &state.mwds[spec_ind][a], lengths[a]-1,diff);
		}
	}
}

/*
 * Converts explicit state representation to compact representation for use 
 * when dumping detailed system contents.
 * Only generates leaves of MWD trees.
 */
void compressState(void) {
	for (int i=0; i < NO_OF_MOLSPECS; i++) {
		int arms = state.arms[i];
		printf("compressing %lld molecules of %s \n", state.ms_cnts[i], name(i));
		if (i >= MAXSIMPLE) {
			// find maximum chain length to determine space required
			chainLen maxLen = 0;
			for (int j=0; j<state.ms_cnts[i]; j++) {
				chainLen *lengths 	= &state.expMols[i].mols[arms*j];
				chainLen totalLen = 0;

				for (int a=0; a<arms; a++) {
					totalLen += lengths[a];
				}
				if (totalLen > maxLen)
					maxLen = totalLen;
			}
			
			// increase allocated memory to store system, if required
			while (maxLen > state.mwds[i][0].maxEntries) {
				printf("%s: mwd tree for %s doubled\n", __FUNCTION__, name(i));
				state.mwds[i][0].maxEntries *= 2;
			}
			int bytes = (2*state.mwds[i][0].maxEntries-1)*sizeof(pcount);
			free(state.mwds[i][0].mwd_tree);
			state.mwds[i][0].mwd_tree = malloc(bytes);
			memset(state.mwds[i][0].mwd_tree,0,bytes);

			// build up MWDs
			for (int j=0; j<state.ms_cnts[i]; j++) {
				chainLen *lengths 	= &state.expMols[i].mols[arms*j];
				chainLen totalLen = 0;

				for (int a=0; a<arms; a++) {
					totalLen += lengths[a];
				}

				IFVERBOSE printf("totalLen = %d maxEntries = %d\n",totalLen,state.mwds[i][0].maxEntries);
				int leavesOffset = state.mwds[i][0].maxEntries-2;
				state.mwds[i][0].mwd_tree[leavesOffset + totalLen]++;
			}
		}
	}
}

/*
 * Converts compact representation to explicit representation.
 */
void decompressState(void) {
	printf("decompressStateRunning... ");
	for (int i=0; i < NO_OF_MOLSPECS; i++) {
		memset(state.expMols[i].mols, 0, sizeof(chainLen)*state.expMols[i].maxMolecules*state.arms[i]);
//		pcount oldMsCnts = state.ms_cnts[i];
		if (i >= MAXPOLY) { // complex species

			break;
			// complex species get slightly more complex treatment since in
			// the compact form, they are represented as pairs of probability
			// distributions.
//			chainLen lens[MAX_ARMS];
//			int dummy;
//			for (int x=0; x < oldMsCnts; x++) {
//				pickRndMolecule(i, lens, &dummy);
//				chainLen *dest = state.expMols[i].mols + x;
//				memcpy(dest, lens, sizeof(chainLen)*state.arms[i]);
//			}
		}
		else if (i >= MAXSIMPLE) { // poly species
			int pos = 0;
			chainLen offset = state.mwds[i][0].maxEntries;
			for (chainLen len=1; len<=offset; len++) {
				while (state.mwds[i][0].mwd_tree[offset-1+len] >=0) {
					state.expMols[i].mols[pos++] = len;
				}
			}
		}
	}
	printf("finished!\n");
}

/* For O(1) algorithm, uses large amounts of memory */
INLINE void adjustMolCnt_order1(int spec_ind, chainLen *lengths, int arms, int diff) {
	pcount molecules = state.ms_cnts[spec_ind];

	if (molecules > state.expMols[spec_ind].maxMolecules) {
		state.expMols[spec_ind].maxMolecules *= 1.5;
		printf("Explicit sys state expanded for species %s from %lld to %lld\n", name(spec_ind), (molecules-1), state.expMols[spec_ind].maxMolecules);
		chainLen *old = state.expMols[spec_ind].mols;
		state.expMols[spec_ind].mols = malloc(state.expMols[spec_ind].maxMolecules*arms*sizeof(chainLen));
		memcpy(state.expMols[spec_ind].mols, old, sizeof(chainLen)*(molecules-1)*arms);
		free(old);
	}
	chainLen *mempos = state.expMols[spec_ind].mols + (molecules-1)*arms;
	for (int i=0; i<arms; i++) {
		mempos[i] = lengths[i];
	}
}

INLINE pcount monomerCount() {
	pcount cnt = 0;
	for (int i = 0; i < MAXMONOMER; i++) {
		cnt += state.ms_cnts[i];
	}
	return cnt;
}

INLINE float conversion(void) {
	float conversion = ((float)state.localMonomerParticles*state.scaleFactor - (float)state.currentMonomerMolecules) / ((float)state.localMonomerParticles*state.scaleFactor);
	return conversion;
}

/*
 *  Returns number of molecules of type spec_ind and length 
 *  (leave_ind +1)
 */
pcount getMolCnt(int spec_ind, int arm, int leave_ind) {
	int ix = leave_ind+state.mwds[spec_ind][arm].maxEntries-1;
	pcount cnt = state.mwds[spec_ind][arm].mwd_tree[ix];

	return (cnt);
}

void dumpReactProbTree() {

	int nextLev = 2;
	
	for (int j=0; j<2*REACT_PROB_TREE_LEAVES-1; j++) {
		printf (" %.40lf",state.reactProbTree[j]);
		if (j == nextLev-2) {
			printf("\n");
			nextLev *= 2;
		}
	}
	printf("\n\n");
}

#if defined(_MSC_VER)
__forceinline void updateTree(int prevReact);
#elif defined(__GNUC__)
inline void updateTree (int prevReact) __attribute__((always_inline));
#endif

INLINE void updateTree(int prevReact) {
	RATES_UPDATE_BODY
#ifdef SIMULATEHEATING
	REACTION_PROBABILITY_TREE_INIT // Use of TREE_UPDATE_BODY isn't sufficient, we need to update the entire tree
#else
	switch (prevReact)
		TREE_UPDATE_BODY
#endif
}

INLINE double toConc(long long ps) {
	return (ps / AVOGADRO / state.volume);
}


/* appends s2 to s1.
* WARNING: unsafe, assumes there is sufficient room for s2 at the end of s1
*/
void strAppend(char *s1, const char *s2) {
	int len1 = strlen(s1);

	if (len1 + strlen(s2) + 1 > MAX_FILENAME_LEN) {
		printf("HARD LIMIT REACHED: MAX_FILENAME_LEN\n");
		printf("%s and %s \n", s1, s2);
		exit(EXIT_FAILURE);
	}
	strcpy(s1 + len1, s2);
}

/* Picks a random reaction by scanning over the rates. 
 *
 * Since the old scan results are still present, it is 
 * not necessary to start recomputing the values from 
 * index 0. It would be sufficient to start with the
 * lowest index of all molecules involved in previous reaction.
 * This would probably pay off for systems with a higher number
 * of reactions.
 *
 * From CodeGen polysim.c.
 */

#if defined(_MSC_VER)
__forceinline int pickRndReaction();
#elif defined(__GNUC__)
inline int pickRndReaction() __attribute__((always_inline));
#endif

/*
* Writes state to files. i.e. for each dist at each time <distname>.<time>.dat in format:
* chain_length\tconc\n
* and <speciesname>.profile.dat for a conc/time profile.
*/
void file_write_state(int mode) {
	int i, j, offset, length;
	char timeStr[MAX_FILENAME_LEN];
	sprintf(timeStr, "%d", (int)roundf(state.time));

	// Initialize writing of concentrations
	char concfname[MAX_FILENAME_LEN] = "\0";
	strAppend(concfname, "concentrations");
	strAppend(concfname, ".csv");
	FILE *conc = fopen(concfname, "a");

	// Check if file could be opened
	if (conc == NULL) {
		printf("Could not write to file %s\n", concfname);
		exit(EXIT_FAILURE);
	}

	// Initialize writing of rates
	char ratesfname[MAX_FILENAME_LEN] = "\0";
	strAppend(ratesfname, "rates");
	strAppend(ratesfname, ".csv");
	FILE *rates = fopen(ratesfname, "a");

	// Check if file could be opened
	if (rates == NULL) {
		printf("Could not write to file %s\n", ratesfname);
		exit(EXIT_FAILURE);
	}

	if (mode == START) {
		// Write headers for concentrations
		fprintf(conc, "Simulation time (min);Conversion");
#ifdef SIMULATEHEATING
		fprintf(conc, ";Temperature (K)");
#endif
		for (i = 0; i < NO_OF_MOLSPECS; i++) {
			fprintf(conc, ";%s (umol/L)", name(i));
		}
#ifdef CALCMOMENTSOFDIST
		fprintf(conc, ";Zeroth moment of dist (umol/L);First moment of dist (umol/L);Second moment of dist (umol/L);Number-average chain length;Weight-average chain length;Polydispersity index");
#endif
		fprintf(conc, "\n");
		
		// Write headers for rates
		fprintf(rates, "Simulation time (min);Conversion");
#ifdef SIMULATEHEATING
		fprintf(rates, ";Temperature (K)");
#endif
		for (i = 0; i < NO_OF_REACTIONS; i++) {
			fprintf(rates, ";%s (1/s)", rname(i));
		}
		fprintf(rates, "\n");

	}
	else if (mode == PROFILES) {

		// Write concentrations
		fprintf(conc, "%f;%f", (state.time / 60.0), state.conversion);
#ifdef SIMULATEHEATING
		fprintf(conc, ";%.2f", state.temp);
#endif
		for (i = 0; i < NO_OF_MOLSPECS; i++) {
			fprintf(conc, ";%e", 1e6*toConc(state.ms_cnts[i]));
		}
#ifdef CALCMOMENTSOFDIST
		fprintf(conc, ";%e", 1e6*toConc(state.zerothMoment));
		fprintf(conc, ";%e", 1e6*toConc(state.firstMoment));
		fprintf(conc, ";%e", 1e6*toConc(state.secondMoment));
		fprintf(conc, ";%e", ((double)state.firstMoment / (double)state.zerothMoment));
		fprintf(conc, ";%e", ((double)state.secondMoment / (double)state.firstMoment));
		fprintf(conc, ";%e", ((double)state.secondMoment * (double)state.zerothMoment / ((double)state.firstMoment * (double)state.firstMoment)));
#endif
		fprintf(conc, "\n");

		// Write rates
		fprintf(rates, "%f;%f", (state.time / 60.0), state.conversion);
#ifdef SIMULATEHEATING
		fprintf(rates, ";%.2f", state.temp);
#endif
		for (i = 0; i < NO_OF_REACTIONS; i++) {
			fprintf(rates, ";%e", state.reactions[i].rc);
		}
		fprintf(rates, "\n");
	}
	else if (mode == FULL) {
	// Write MWD of each poly/complex species
		for (i = 0; i < NO_OF_MOLSPECS; i++) {
			if (isPolyMol(i) && mode == FULL) {
				char distfname[MAX_FILENAME_LEN] = "\0";
				strAppend(distfname, name(i));
				strAppend(distfname, "-");
				strAppend(distfname, timeStr);
				strAppend(distfname, ".csv");
				FILE *dist = fopen(distfname, "w");

				// Check if file could be opened
				if (dist == NULL) {
					printf("Could not write to file %s\n", distfname);
					exit(EXIT_FAILURE);
				}

				fprintf(dist, "Chain length;Concentration (umol/L)\n");
				offset = state.mwds[i][0].maxEntries - 1;
				length = 1;
				for (j = offset; j < 2 * offset + 1; j++) {
					if (state.mwds[i][0].mwd_tree[j] > 0) {
						fprintf(dist, "%d;%e\n", length, 1e6*toConc(state.mwds[i][0].mwd_tree[j]));
					}
					length++;
				}
				fclose(dist);
			}
		}
	}
	fclose(conc);
	fclose(rates);
}

INLINE int pickRndReaction() {
    react_prob     rate;
    int             i;
    react_prob          rnd, origRnd;
	static int prevReact;

	// Choose reaction based on probability
	updateTree(prevReact);
	rate = state.reactProbTree[0];
    rnd = randomProb(2);

	rnd *= rate;
	origRnd = rnd;

	int curr = 1;
	react_prob curr_prob;
  	while (curr < REACT_PROB_TREE_LEAVES) {
		curr = curr << 1;
		curr_prob = state.reactProbTree[curr - 1];
		if (rnd >= curr_prob) {
			curr++;
			rnd -= curr_prob;
		}
  	}
	i = curr - REACT_PROB_TREE_LEAVES;
	prevReact = i;

	if (rate <= 0) {

		RANK compressState();
		RANK file_write_state(FULL);
		if (numprocs == 1) {
			printf("Fatal error, rate = %lf, dumping system and exiting\n", rate);
			print_state();
			if (rate < 0) {
				dumpReactProbTree();
			}
			exit(EXIT_FAILURE);
		}
		else {
			return(-1);
		}
	}

#if DEBUG
    if (i >= NO_OF_REACTIONS) {
        printf("pickRndReaction: ran out of reactions, i = %d no reactions = %d\n",i,NO_OF_REACTIONS);

		dumpReactProbTree();
		REACTION_PROBABILITY_TREE_INIT
		dumpReactProbTree();

		printf("\n original rnd no was %.40lf\n",origRnd);
        exit(EXIT_FAILURE);
    }
#endif

#ifdef TRACE
    printf("reaction %d: ", i);
    print_reaction(i);
#endif

    return i;
}

/* Reads pairs of numbers from data file in format "%f %f" and stores them
 * in struct p.
 */
void readDataFile(const char *fname, TimesVals *p) {
	FILE *f;
	if ((f = fopen(fname,"r")) == NULL) {
		printf("Couldn't open data file %s\n",fname);
		exit (1);
	}

	p->ix = 0;
	p->maxIx = 0;

	float time,value;
	while (fscanf(f,"%f %f\n", &time, &value) != EOF 
			&& p->maxIx < MAX_DATA_FILE)
	{
		p->ts[p->maxIx] = time;
		p->xs[p->maxIx++] = value;
	}

	if (p->maxIx >= MAX_DATA_FILE) {
		printf("Increase MAX_DATA_FILE\n");
		exit(EXIT_FAILURE);
	}
	fclose(f);

	printf("%s: %d (time,x) pairs read from file %s\n",__FUNCTION__,p->maxIx,fname);

}

INLINE int MPI_Recv_int_wrap() {
	int i;
	const int howmany = 1;
	MPI_Bcast(&i, howmany, MPI_INT, 0, MPI_COMM_WORLD);
	return i;
}

void setupParallelData() {
	if (myid == 0) {
		int msg;
		// send SETUP_CONVDATA
		msg = SETUP_CONVDATA;
		MPI_Bcast (&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);

		msg = state.timeCalcData.maxIx;
		MPI_Bcast (&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast (state.timeCalcData.ts, state.timeCalcData.maxIx, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast (state.timeCalcData.xs, state.timeCalcData.maxIx, MPI_FLOAT, 0, MPI_COMM_WORLD);
		printf("Master: broadcasted moments data\n");

		// send SETUP_SYSTEMSCALES
#ifdef SCALING
		msg = SETUP_SYSTEMSCALES;
		MPI_Bcast (&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);

		msg = state.sysScaleData.maxIx;
		MPI_Bcast (&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast (state.sysScaleData.ts, state.sysScaleData.maxIx, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast (state.sysScaleData.xs, state.sysScaleData.maxIx, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
		
		// Send SETUP_END
		msg = SETUP_END;
		MPI_Bcast (&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	else {
		printf("Slave receiving...\n");
		for(int task = MPI_Recv_int_wrap(); task != SETUP_END; task = MPI_Recv_int_wrap()) {
			if (task == SETUP_CONVDATA) {
				int howmany;
			    MPI_Bcast(&howmany, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(state.timeCalcData.ts, howmany, MPI_FLOAT, 0, MPI_COMM_WORLD);
				MPI_Bcast(state.timeCalcData.xs, howmany, MPI_FLOAT, 0, MPI_COMM_WORLD);
				state.timeCalcData.maxIx = howmany;
				printf("%d moments data points received from Master\n",state.timeCalcData.maxIx);
			}
			else if (task == SETUP_SYSTEMSCALES) {
				int howmany;
			    MPI_Bcast(&howmany, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(state.sysScaleData.ts, howmany, MPI_FLOAT, 0, MPI_COMM_WORLD);
				MPI_Bcast(state.sysScaleData.xs, howmany, MPI_FLOAT, 0, MPI_COMM_WORLD);
				state.sysScaleData.maxIx = howmany;
				printf("%d system scales received from Master\n",state.sysScaleData.maxIx);
			}
			else {
				printf ("Error: bad setup task\n");
				exit(EXIT_FAILURE);
			}
		}
	}
}

/* Initialise the state variable.
 * 
 * From CodeGen polysim.c.
 */
void initSysState(int seed) {
    int i;

	dsfmt_gv_init_gen_rand(seed);

#ifdef SCALING
	RANK readDataFile("system.scaling",&state.sysScaleData);
#endif
	state.scaleFactor = 1.0;

	if (numprocs > 1)
		setupParallelData();

    /* time, rate, particles, events, conversion, temperatures, etc. */
    state.time = 0.0;
	state.nextSynchTime = state.synchTime;
	state.events = 0;
	state.noMoreReactions = 0;
	state.nextSynchEvents = (int)((float)state.synchEvents/(float)numprocs);
	state.localMonomerParticles = takeSome(GLOBAL_MONOMER_PARTICLES);
	state.conversion = 0;
	state.basetemp = BASETEMP;
	state.temp = STARTTEMP;
	state.deltatemp = state.temp - state.basetemp;
	
    /* reactions. used for temporary memory in pickRndReactions */
    state.scan_scratch = calloc(NO_OF_REACTIONS, sizeof(probability));

#ifdef CALCMOMENTSOFDIST
	/* moments of distribution */
	state.zerothMoment = 1;
	state.firstMoment = 1;
	state.secondMoment = 1;
#endif
#ifdef CALCFREEVOLUME
	/* free volume */
	state.freeVolumeFraction = (VF0 + ALPHA_P * (state.temp - TG_P)) * state.conversion + (VF0 + ALPHA_M * (state.temp - TG_M)) * (1 - state.conversion);
#endif
    /*
     * MWDS initialise those ms_cnts not set to 0. rely on calloc for the
     * rest
     */
    MWD_INITS
	
	for (i = 0; i < NO_OF_MOLSPECS; i++) {
		state.ms_cnts[i] = takeSome(state.ms_cnts[i]);
	}

    /* initialise mwds */
    /* enter all mols into mwds */
	printf("startsize = %d\n",startsize);
    for (i = 0; i < NO_OF_MOLSPECS; i++) {
		for (int a = 0; a < MAX_ARMS; a++) {
	        state.mwds[i][a].maxEntries = startsize;
	        state.mwds[i][a].mwd_tree = calloc(2 * startsize - 1, sizeof(pcount));
		}
    }

    //state.reactions = malloc(sizeof(reaction) * NO_OF_REACTIONS);

	REACTIONS_INIT

	// Adjust rate coefficients for bimolecular reactions.
	// Needed because particles might be divided between > 1 node.
	for (i = 0; i < NO_OF_REACTIONS; i++) {
		if (state.reactions[i].arg_ms2 != NO_MOL) {
			state.reactions[i].rc *= ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles);
		}
	}

	state.initialMonomerMolecules = monomerCount();
	state.currentMonomerMolecules = state.initialMonomerMolecules;
	state.volume = (state.localMonomerParticles/AVOGADRO)/CONCENTRATION;

	int tmp[NO_OF_MOLSPECS] = ARMS_INIT;
	memcpy (state.arms, tmp, NO_OF_MOLSPECS*sizeof(int));
	memset (state.reactProbTree, 0, (2*REACT_PROB_TREE_LEAVES-1)*sizeof(probability));
	REACTION_PROBABILITY_TREE_INIT

	
#ifdef EXPLICIT_SYSTEM_STATE
	state.expMols = malloc(sizeof(MoleculeList)*NO_OF_MOLSPECS);
	for (int i=0; i<NO_OF_MOLSPECS; i++) {
		int initialSize = 10485760;
		if (i >= MAXSIMPLE) {
			state.expMols[i].mols = malloc(initialSize*sizeof(chainLen)*state.arms[i]);
			state.expMols[i].maxMolecules = initialSize;
		}
	}
#endif

}

//#if defined(_MSC_VER)
//__forceinline int pickRndChainLen(pcount *mwd_tree, int size);
//#elif defined(__GNUC__)
//inline int pickRndChainLen(pcount *mwd_tree, int size) __attribute__((always_inline));
//#endif

INLINE int pickRndChainLen(pcount *mwd_tree, int size) {

	probability     prob, curr_prob;

	int             length;
	int             curr = 1;

	prob =  randomProb(3);  
	prob *= mwd_tree[0];

	mwd_tree[0]--;
	while (curr < size) {
		curr = curr << 1;
		curr_prob = (probability)mwd_tree[curr - 1];
		if (prob >= curr_prob) {
			curr++;
			prob -= curr_prob;
		}
#ifdef DEBUG
		if (mwd_tree[curr - 1] < 0) {
			return (-1);
		}
#endif
		mwd_tree[curr - 1]--;
	}

	length = curr - size + 1;
	return (length);

}

/*
 *  Pick a random molecule which matches the spec and delete it from
 *  system. Returns: length of molecule picked.
 * From CodeGen polysim.c.
 */
int pickRndMolecule(int spec_index, chainLen *lens, int *arms) {

#if defined(DEBUG)
    //assert(state.ms_cnts[spec_index] > 0);
	if (state.ms_cnts[spec_index] <= 0) {
		printf("Trying to pick %s when there are none\n",name(spec_index));
		print_state();
		dumpReactProbTree();
		return -1;
	}
#endif

    /* decrement global and spec local count */
//    state.no_of_mols--;
	if (state.ms_cnts[spec_index] > 0)
	    state.ms_cnts[spec_index]--;

    /* simple molecule, we're done */
    if (spec_index < MAXSIMPLE) {
        return 0;
    }
	else if (spec_index < MAXPOLY) {
		/* Poly species: pick one chain lengths according to MWD */
		*arms = 1;
		pcount *mwd_tree = state.mwds[spec_index][0].mwd_tree;
		int size         = state.mwds[spec_index][0].maxEntries;
		lens[0] = pickRndChainLen(mwd_tree,size);
#ifdef DEBUG
		if (lens[0] < 0) {
			printf("error: chain length < 0 for species %s(arm=0)\n", name(spec_index));
			print_state();
			dumpReactProbTree();
			abort();
		}
#endif

	}
	else {

		/* Complex species: pick several chain lengths according to MWD */
		*arms = state.arms[spec_index];
		for (int a=0; a<*arms; a++) {
			pcount *mwd_tree = state.mwds[spec_index][a].mwd_tree;
			int size         = state.mwds[spec_index][a].maxEntries;
			lens[a] = pickRndChainLen(mwd_tree,size);
#ifdef DEBUG
			if (lens[a] < 0) {
				printf("error: chain length < 0 for species %s(arm=%d)\n", name(spec_index), a);
				print_state();
				dumpReactProbTree();
				abort();
			}
#endif
		}
	}

	return 0;
}

int pickRndMolecule_order1(int spec_index, chainLen *lens, int *arms) {

#if defined(DEBUG)
    //assert(state.ms_cnts[spec_index] > 0);
	if (state.ms_cnts[spec_index] <= 0) {
		printf("Trying to pick %s when there are none\n",name(spec_index));
		print_state();
		dumpReactProbTree();
		return -1;
	}
#endif

    /* simple molecule, we're done */
    if (spec_index < MAXSIMPLE) {
    	state.ms_cnts[spec_index]--;
        return 0;
    }
	else {
		*arms = state.arms[spec_index];
		int lastIx = state.ms_cnts[spec_index]-1;
        // int rndIx = (int)round((double)randomProb(3)*(double)lastIx);
		/* Fast rounding by adding 0.5, then casting to int
		   In theory, this may yield incorrect results in cases where 'x' and/or 'x + 0.5' is not exactly representable and 
		   x is very close to a half-integer. In practice, this will virtually never happen. */
		int rndIx = (int)((double)randomProb(3)*(double)lastIx + 0.5);
		chainLen *pickedPos = &state.expMols[spec_index].mols[rndIx*(*arms)];
		int memSize = (*arms)*sizeof(chainLen);
		memcpy(lens, pickedPos, memSize);

		if (rndIx != lastIx) {
			memcpy(pickedPos, &state.expMols[spec_index].mols[lastIx*(*arms)], memSize);
		}

	    /* decrement global and spec local count */
	    state.ms_cnts[spec_index]--;
		return 0;
	}
}

INLINE const char * name(int index) {
    switch (index) {
        MOLECULENAMES
    }
    return 0;
}


INLINE const char * rname(int index) {
	switch (index) {
		REACTIONNAMES
	}
	return 0;
}

/* Scales number of particles and (rate coefficients accordingly) by
 * factor.
 */
void scaleSystem(float factor) {

	if (factor < 1.0F) {
		printf("Scaling factors < 1 not yet dealt with\n");
		abort();
	}
	
	print_state();
	state.scaleFactor = factor;

	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		if (isSimpleSpec(i)) {
			state.ms_cnts[i] *= factor;
			continue;
		}

		for (int a=0; a < state.arms[i]; a++) {
			int j;
			for (j=0; j < state.mwds[i][a].maxEntries-1; j++) {
				state.mwds[i][a].mwd_tree[j] *= factor;
			}

			pcount total = 0;
			for (j=state.mwds[i][a].maxEntries-1; j < state.mwds[i][a].maxEntries*2-1; j++) {
				state.mwds[i][a].mwd_tree[j] *= factor;
				total += state.mwds[i][a].mwd_tree[j];
			}
//			printf("System scale: %s changed to %lld particles\n",name(i),total);
			state.ms_cnts[i] = total; // total molecules are counted from tree to avoid rounding errors
		}
	}

	// Scale rate coefficients.
	for (int i = 0; i < NO_OF_REACTIONS; i++) {
		if (state.reactions[i].arg_ms2 != NO_MOL) { // only if it's not unimolecular
			state.reactions[i].rc /= factor;
		}
	}

	state.initialMonomerMolecules *= factor;
	state.currentMonomerMolecules *= factor;
	state.volume *= factor;

	printf("System scaled by factor of %f\n",factor);
	print_state();

}

/* Print reactants and results of a reaction
 *   (we could print the name of the reaction as well, would need to
 *   generate code for this -- see name (..))
 */
void print_reaction(int reaction_index) {
    printf("%s ", name(state.reactions[reaction_index].arg_ms1));
    if (state.reactions[reaction_index].arg_ms2 < NO_MOL) {
        printf("+ %s ", name(state.reactions[reaction_index].arg_ms2));
    }
    printf("\t--> %s ", name(state.reactions[reaction_index].res_ms1));
    if (state.reactions[reaction_index].res_ms2 < NO_MOL) {
        printf("+ %s ", name(state.reactions[reaction_index].res_ms2));
    }
    printf("\n");
}

/* Check if tree is still consistent:
 *  - does not contain negative values in any node
 *  - sum of leave value equals node value for all nodes,
 *    all leaves
 */
pcount check_tree(int no_of_leaves,
		int levelSize,
        int offset,
        pcount * mwd,
	    int spec_ind)
{
    pcount     left, right;

    if (levelSize == no_of_leaves) {
        if (mwd[offset] >= 0) {
            return (mwd[offset]);
        } else {
            printf("check_tree: negative value in tree for species %s\n", name(spec_ind));
            print_state();
            exit(EXIT_FAILURE);
        }
    } else {
		int ixInLevel = 2*offset+1;
        left  = check_tree(no_of_leaves, 
							levelSize*2,
                  			ixInLevel,
							mwd,
							spec_ind);
        right = check_tree(no_of_leaves,
							levelSize*2,
                   			ixInLevel+1, 
							mwd,
							spec_ind);
        if (mwd[offset] == left + right) {
			return (mwd[offset]);
        } else {
			printf("check_tree: bad tree for species %s\n", name(spec_ind));
			printf("left = %lld right = %lld sum = %lld offset = %d",left,right,mwd[offset],offset);
            print_state();
            exit(EXIT_FAILURE);
        }
    }
}

/* React body. Optimized for performing as many iterations in as little time as possible.
 * For SYNCH_MODE == 0, a 'do while' loop is the fastest
 * For SYNCH_MODE == 1, a 'for' loop is the fastest
 * (Tested with MSVC, should be tested also with GCC)
 */
void react(void) {

	int             reactionIndex;

	int             react1_ind, react2_ind;
	chainLen        react1_lens[MAX_ARMS], react2_lens[MAX_ARMS];
	int				react1_arms = 1, react2_arms = 1, prod1_arms = 1, prod2_arms = 1;

	int				prod1_ind, prod2_ind;
	chainLen		prod1_lens[MAX_ARMS], prod2_lens[MAX_ARMS];

	int				no_of_res;

	int				prm1, prm2;
	
	int chainLenLimit = (int)(pow(2, (8 * sizeof(chainLen))) - 0.5);

#if SYNCH_MODE == 0
	do {
#elif SYNCH_MODE == 1
	for (; state.events < state.nextSynchEvents; state.events++) {
#endif

		/* basic case */
		react1_arms = 1;
		react2_arms = 1;
		prod1_arms = 1;
		prod2_arms = 1;
		react2_ind = -1;
		prod2_ind = -1;
		no_of_res = 1; 
		prm1 = 0;
		prm2 = 0;

		reactionIndex = pickRndReaction();

		// No events are possible, system does nothing
		if (reactionIndex == -1) {
			state.time = state.nextSynchTime;
			state.events = state.nextSynchEvents;
			state.noMoreReactions = 1;
			return;
		}

		react1_ind = reactToSpecInd1(state.reactions[reactionIndex]);

#ifdef EXPLICIT_SYSTEM_STATE
		prm1 = pickRndMolecule_order1(react1_ind, react1_lens, &react1_arms);
#else
		prm1 = pickRndMolecule(react1_ind, react1_lens, &react1_arms);
#endif

		if (SPECTYPE(state.reactions[reactionIndex].arg_ms2) != NO_MOL) {
			react2_ind = reactToSpecInd2(state.reactions[reactionIndex]);
#ifdef EXPLICIT_SYSTEM_STATE
			prm2 = pickRndMolecule_order1(react2_ind, react2_lens, &react2_arms);
#else
			prm2 = pickRndMolecule(react2_ind, react2_lens, &react2_arms);
#endif

		}

#ifdef DEBUG
		if (prm1 != 0 || prm2 != 0) {
			printf("Problem reaction %d\n", reactionIndex);
			abort();
		}
#endif

		/*
		 * Given a reaction index, the lengths l1 and l2  of the reactants (0
		 * in case of Simple), determine type and length(s) of resulting
		 * molecule(s)
		 */
		switch (reactionIndex)
			DO_REACT_BODY

#ifdef DEBUG
		assert(prod1_ind != -1);
#endif

		state.ms_cnts[prod1_ind]++;
		if (prod1_ind >= MAXSIMPLE) {
#ifdef EXPLICIT_SYSTEM_STATE
			adjustMolCnt_order1(prod1_ind, prod1_lens, prod1_arms, 1);
#else
			adjustMolCnt(prod1_ind, prod1_lens, prod1_arms, 1);
#endif

#ifdef DEBUG
			assert(prod1_arms == state.arms[prod1_ind]);
#endif
		}
		if (no_of_res == 2) {
			state.ms_cnts[prod2_ind]++;
			if (prod2_ind >= MAXSIMPLE) {
#ifdef EXPLICIT_SYSTEM_STATE
				adjustMolCnt_order1(prod2_ind, prod2_lens, prod2_arms, 1);
#else
				adjustMolCnt(prod2_ind, prod2_lens, prod2_arms, 1);
#endif
			}
		}

		// Integer overflow protection
		if ((prod1_lens[0] == chainLenLimit) ||
			(prod2_lens[0] == chainLenLimit)) {
			printf("\nError: a species on node %d has reached its maximum chain length.\nWriting data and exiting...\n\n", myid);
			RANK compressState();
			RANK file_write_state(FULL);
			printf("\nNode %d exited with errors\n", myid);
			exit(EXIT_FAILURE);
		}

#ifdef CALCMOMENTSOFDIST
		// Adjust moments for disappearance reactant 1
		if (react1_ind >= MAXSIMPLE) {
			state.zerothMoment -= 1;
			if (react1_ind < MAXPOLY) {
				state.firstMoment -= react1_lens[0];
				state.secondMoment -= react1_lens[0] * react1_lens[0];
			}
			else {
				int totalLength = 0;
				for (int arm = 0; arm < react1_arms; arm++) {
					totalLength += react1_lens[arm];
				}
				state.firstMoment -= totalLength;
				state.secondMoment -= totalLength * totalLength;
			}
			// Adjust moments for disappearance reactant 2
		}
		if (react2_ind >= MAXSIMPLE) {
			state.zerothMoment -= 1;
			if (react2_ind < MAXPOLY) {
				state.firstMoment -= react2_lens[0];
				state.secondMoment -= react2_lens[0] * react2_lens[0];
			}
			else {
				int totalLength = 0;
				for (int arm = 0; arm < react2_arms; arm++) {
					totalLength += react2_lens[arm];
				}
				state.firstMoment -= totalLength;
				state.secondMoment -= totalLength * totalLength;
			}
		}
		// Adjust moments for appearance product  1
		if (prod1_ind >= MAXSIMPLE) {
			state.zerothMoment += 1;
			if (prod1_ind < MAXPOLY) {
				state.firstMoment += prod1_lens[0];
				state.secondMoment += prod1_lens[0] * prod1_lens[0];
			}
			else {
				int totalLength = 0;
				for (int arm = 0; arm < prod1_arms; arm++) {
					totalLength += prod1_lens[arm];
					state.firstMoment += totalLength;
					state.secondMoment += totalLength * totalLength;
				}
			}
		// Adjust moments for appearance product 2
		}
		if (prod2_ind >= MAXSIMPLE) {
			state.zerothMoment += 1;
			if (prod2_ind < MAXPOLY) {
				state.firstMoment += prod2_lens[0];
				state.secondMoment += prod2_lens[0] * prod2_lens[0];
			}
			else {
				int totalLength = 0;
				for (int arm = 0; arm < prod2_arms; arm++) {
					totalLength += prod2_lens[arm];
					state.firstMoment += totalLength;
					state.secondMoment += totalLength * totalLength;
				}
			}
		}
#endif

#ifdef RECALCCONVERSION
		// Blindly recalculating the conversion is cheaper than first checking whether the conversion needs to be recalculated
		// pcount prevMonomerCount = state.currentMonomerMolecules;
		state.currentMonomerMolecules = monomerCount();
		//if (prevMonomerCount != state.currentMonomerMolecules) {
		state.conversion = conversion();
		//}
#endif

#ifdef SIMULATEHEATING
		state.deltatemp += (state.reactions[reactionIndex].energy / state.volume);
		state.temp = state.basetemp + state.deltatemp;
#endif

#if defined(TRACE)
		if (monomerAudit("Trace") == -1)
			print_state();
		//	dumpReactProbTree();
		printf("--------------------------\n");
#endif

#ifdef MICRO_CALC_TIME
		float           rndtime;
		probability     rate;
		ptime			deltatime;

		rndtime = randomProb(2);
		rate = (float)state.reactProbTree[0];

		deltatime = (-logf(rndtime)) / rate;

#ifdef COOLINGRATE
		if (state.deltatemp != 0.0) { // Don't cool when not possible
			state.deltatemp *= exp(-1 * COOLINGRATE * deltatime);
			state.temp = state.basetemp + state.deltatemp;
		}
#endif

		state.time += deltatime;
#endif

#ifdef CALCFREEVOLUME
		state.freeVolumeFraction = (VF0 + ALPHA_P * (state.temp - TG_P)) * state.conversion + (VF0 + ALPHA_M * (state.temp - TG_M)) * (1 - state.conversion);
#endif

#if defined(DEBUG)
		//	monomerAudit(MONOMER_INDEX);
		assert(rate >= 0);
#endif

#if SYNCH_MODE == 0
		state.events++;
#endif
	}
#if SYNCH_MODE == 0
	while (state.time < state.nextSynchTime);
#endif
}

// Prints the mwds and molecule counts
 void print_state() {
    int             i, j, offset;

    for (i = 0; i < NO_OF_MOLSPECS; i++) {
        if (state.ms_cnts[i] > 0) {
            printf("\n\n# %s (%lld):\n", name(i), state.ms_cnts[i]);
            if (isSimpleSpec(i)) {
                printf("%lld\n", state.ms_cnts[i]);
            } else {
				pcount cnts[MAX_ARMS];
				offset = state.mwds[i][0].maxEntries;
				for (j = 0; j < offset; j++) {
					int total = 0;
					for (int a=0; a<state.arms[i]; a++) {
						cnts[a] = getMolCnt(i, a, j);
						total += cnts[a];
					}
					if (total > 0) {
						printf("%d",j+1);
						for (int a=0; a<state.arms[i]; a++) {
							printf("\t%lld", cnts[a]);
						}
						printf("\n");
					}
				}
            }
        }
    }
}

void print_state_summary(void) {

	int nameLens[NO_OF_MOLSPECS];
	int maxNumLen = strlen("0.0e+00");
	int sum = 0;

	for (int i=0; i < NO_OF_MOLSPECS; i++) {
		nameLens[i] = max_value(strlen(name(i)), maxNumLen) + 1;
		sum += nameLens[i];
	}

	for (int j=0; j<sum; j++)
		printf("-");
	printf("\n");

	printf("Time (Wall, Poly) = %lf %f\n", readTimerSec(&state.wallTime), state.time);

#ifdef SIMULATEHEATING
	printf("Temperature = %.2f K\n", (roundf(state.temp * 100) / 100));
#endif

	for (int i=0; i < NO_OF_MOLSPECS; i++) {
		printf("%s",name(i));
		for (int j = 0; j < (nameLens[i]-strlen(name(i))); j++)
			printf(" ");
	}
	printf("\n");

	for (int i=0; i < NO_OF_MOLSPECS; i++) {
		char tmp[MAX_FILENAME_LEN];
		sprintf(tmp, "%.1E", 1e6*toConc(state.ms_cnts[i]));
		printf("%s",tmp);
		for (int j = 0; j < (nameLens[i]-strlen(tmp)); j++)
			printf(" ");
	}
	printf("\n");

	for (int j=0; j<sum; j++)
		printf("-");
	printf("\n");

}

void print_state_conc(void) {
  int i, j, offset, length;
 
  printf ("System state:\n\ttime: %f node volume = %E L\n",state.time,state.volume);
   
  for (i = 0; i < NO_OF_MOLSPECS; i++) {
    if (state.ms_cnts[i] > 0) {
      printf (" %s %lld = (%.0f umol/L):\n", name(i), state.ms_cnts[i], 1e6*toConc(state.ms_cnts[i]));
      if (isSimpleSpec(i)) {
		printf ("\t %lld = %.0f umol/L\n", state.ms_cnts[i], 1e6*toConc(state.ms_cnts[i]));
      } else {
	offset = state.mwds[i][0].maxEntries - 1;
	length = 1;
	for (j = offset; j < 2*offset+1; j++) {
	  if (state.mwds[i][0].mwd_tree[j] > 0) {
		printf ("%d\t%lld = %.0f umol/L \n", length, state.mwds[i][0].mwd_tree[j], 1e6*toConc(state.mwds[i][0].mwd_tree[j]));
	  }
	  length++;
	}
      }
    }
  }
  
}

void print_kinetic_model(void) {
	printf("\tname\n");
	for (int i=0; i<NO_OF_MOLSPECS; i++) {
        char *typeStr;
			if (isPolyMol(i))
				typeStr = "(i)";
			else if (isComplexMol(i))
				typeStr = "(i,j,...)";
			else
				typeStr = "";

		printf("%d:\t%s%s\n",i,name(i),typeStr);
	}
	printf("\nKinetic model\n");
	for (int i=0; i<NO_OF_REACTIONS; i++) {
		reaction *p = state.reactions+i;
		char *p1 = (char*)((p->res_ms1!=NO_MOL)?name(p->res_ms1):"  ");
		char *p2 = (char*)((p->res_ms2!=NO_MOL)?name(p->res_ms2):"  ");
		char *r1 = (char*)((p->arg_ms1!=NO_MOL)?name(p->arg_ms1):"  ");
		char *r2 = (char*)((p->arg_ms2!=NO_MOL)?name(p->arg_ms2):"  ");
		double k = state.reactions[i].rc;
		printf("%d:\t%s + %s\t-->\t%s + %s\t k = %.5e\n",i,r1,r2,p1,p2,k);
	}

}

pcount getConvertedMonomer() {
    int i, j;
	pcount convertedMonomer = 0;

#ifdef EXPLICIT_SYSTEM_STATE
    for (i = 0; i < NO_OF_MOLSPECS; i++) {
		if (i >= MAXSIMPLE) {
			for (j = 0; j < state.arms[i]*state.ms_cnts[i]; j++) {
				convertedMonomer += state.expMols[i].mols[j];
			}
		}
    }
#else
	int length, offset;
    for (i = 0; i < NO_OF_MOLSPECS; i++) {
		if (!isSimpleMol(i)) {
			for (int a = 0; a < state.arms[i]; a++) {
				length = 0;
				offset = state.mwds[i][a].maxEntries - 1;
				for (j = offset; j < 2*offset+1; j++) {
					length++;
					convertedMonomer += length*state.mwds[i][a].mwd_tree[j];
				}
			}
		}
    }
#endif

	return convertedMonomer;
}

int monomerAudit(const char *str) {

	printf("Auditing on node %d...\n", myid);
    pcount discrepancy = state.currentMonomerMolecules
					   + getConvertedMonomer()
					   - state.initialMonomerMolecules;

    if (discrepancy == 0) {
	    printf("Local monomer audit (%s) passed!!!\n", str);
		return (0);
	}
    printf("Local monomer audit (%s) failed!!! Discrepancy = %lld\n", str, discrepancy);
	return (-1);
}


// ************* begin communication code ********

int stateCommSize = INIT_STATE_COMM_SIZE;

void dumpStateComm(StatePacket *x,char *s) {
	printf("stateComm in %s...\n",s);
	printf("sc.tv_sec = %f sc.tooBig = %d\n",x->time,x->stateTooBig);
	printf("comm|state\n");

	pcount *molCounts = (pcount*)(x+1);
	
	int *maxChainLens = (int*)(molCounts+NO_OF_MOLSPECS);
//	int tmp;

	int maxMolCnts = 0;
	for (int i=0; i < NO_OF_MOLSPECS; i++) {
		maxMolCnts += state.arms[i];
	}
	printf("molCounts:");
	for (int i=0; i < NO_OF_MOLSPECS; i++) {
		printf(" %lld|%lld",molCounts[i],state.ms_cnts[i]);
	}
	printf("\nmaxChainLens");

	for (int i=0; i < maxMolCnts; i++) {
		printf(" %d|%d",maxChainLens[i],state.mwds[i][0].maxEntries);
	}
	pcount *mwd = (pcount*)(maxChainLens+maxMolCnts);
	printf("\nsample chain lengths");
/*	for (int j=0; j < maxMolCnts; j++) {
		for (int i=0; i< maxChainLens[j]; i++) {
			printf(" %lld",mwd[i]);
		}
	}*/
	for (int i=0; i< 30; i++) {
        printf(" %lld",mwd[i]);
    }
	printf("...\n");
	
/*
	pcount *mwdTrees = (pcount*)(mwdTreeMaxEntries+NO_OF_MOLSPECS);
	for (int i=0; i < NO_OF_MOLSPECS; i++) {
		printf("\n%d:");

		for (int j=0; j<TREE_SAMPLE;j++) {
			printf(" %ld|%ld",mwdTrees[j],state.mwds[i].mwd_tree[j]);
			
		}
		tmp = mwdTreeMaxEntries[i];
		mwdTrees += tmp;
    }
*/	
	printf("\n");
}


/* Performs a+b=c for array of doubles
 */
INLINE void dvecAdd(pcount *c, pcount *a,pcount *b,int n) {
	for (int i=0; i<n; i++) {
			c[i] = a[i]+b[i];
	}
}

INLINE int sizeExceeded(StatePacket *sp) {
    return(sp->stateTooBig);
}


void buildTreeFromLeaves(pcount *tree, pcount *leaves, int maxLeaves) {
	//safety check
	if (maxLeaves % 2 != 0) {
		printf ("buildTreesFromLeaves had uneven leaves count\n");
		exit(EXIT_FAILURE);
	}

	// stores cummulative sum of MWD, i.e. the integral
	pcount *sums = malloc ((1+maxLeaves) * sizeof(pcount));
	sums[0] = 0;
	for (int i=1; i < maxLeaves+1; i++) {
		sums[i] = sums[i-1] + leaves[i-1];
	}
	
	int j = 0; // index into new tree as it is built
	for (int subTreeSize = maxLeaves; subTreeSize > 1; subTreeSize /= 2) {
		for (int i=0; i < maxLeaves; i+=subTreeSize) {
			tree[j++] = sums[i+subTreeSize] - sums[i];
		}
	}
	free(sums);
	//safety check
	if (j != maxLeaves-1) {
		printf("buildTreesFromLeaves failed final safety check j = %d maxLeaves = %d\n",j,maxLeaves);
		exit(EXIT_FAILURE);
	}
}

/*
 * populates list of the max chain length of all species. Simple 
 * species have chain length = 0. returns bytes needed to store all
 * MWDs.
 */
int getMaxChainLens(int *max_len_arr_ptr) {
	int bytesNeeded = 0, pos = 0;
	
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		if (isSimpleMol(i)) {
			max_len_arr_ptr[pos] = 0;
			pos++;
		} else {
#ifdef EXPLICIT_SYSTEM_STATE
			if (state.arms[i] == 1) { // for single arm species, stores max chain length as usual
				chainLen maxLen = 0;
				for (int x=0; x<state.expMols[i].maxMolecules; x++) {
					chainLen tmp = state.expMols[i].mols[x];
					if (tmp > maxLen)
						maxLen = tmp;
				}
				max_len_arr_ptr[pos] = maxLen;
				IFVERBOSE printf("max_len_arr_ptr[%s] = %d\n",name(i),maxLen);
				bytesNeeded += max_len_arr_ptr[pos];
				pos++;
			}
			else { // for complex species, 1st entry stores number of molecules, rest are zero
				max_len_arr_ptr[pos] = state.ms_cnts[i];
				IFVERBOSE printf("max_len_arr_ptr[%s] = %d\n",name(i),max_len_arr_ptr[pos]);
				pos++;
				for (int a=1; a<state.arms[i]; a++) {
					max_len_arr_ptr[pos++] = 0;
				}
			}
#else
			int currMax;
			for (int a=0; a<state.arms[i]; a++) {
				int offset = state.mwds[i][a].maxEntries - 1;
				currMax = offset-1;
				for (int j = offset; j < 2*offset+1; j++) {
					if (state.mwds[i][a].mwd_tree[j] > 0) {
						currMax = j;
					}
				}
				max_len_arr_ptr[pos] = currMax-offset+1;
				bytesNeeded += max_len_arr_ptr[pos];

				pos++;
			}
#endif
		}
    }
	bytesNeeded *= sizeof(pcount);
	return bytesNeeded;
}

int stateCommHeaderBytes(void) {
	int maxChainLensSize = 0; // size required to store max len of every species in sys
	for (int s=0; s < NO_OF_MOLSPECS; s++) {
		maxChainLensSize += state.arms[s];
	}

	return (sizeof(StatePacket)
		  + sizeof(pcount)*NO_OF_MOLSPECS // mol counts
          + sizeof(int)*maxChainLensSize);  // max chain lengths
}

int requiredStateCommSize(void) {

	int chainLens = 0;
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		chainLens += state.arms[i];
	}

	int *maxChainLens = malloc(chainLens * sizeof(int));
	int bytesForMWDs = getMaxChainLens(maxChainLens);

	free(maxChainLens);

	return stateCommHeaderBytes() + bytesForMWDs;
}

void stateToComm(StatePacket **outStatePacket, StatePacket **inStatePacket) {

	// Packet format
	// [ StatePacket    | pcount[]  | int[]        | pcount[] | pcount[] | ... | pcount[]]
	// [ outStatePacket | molCounts | maxChainLens | packetMWDs ...                    E ]
	// outStatePacket is statepacket struct
	// molCounts      is number of each species (NO_OF_MOLSPECS elts)
	// maxChainLens   is the max chainlength stored for each dist, or MIN_MWD_SIZE, 
    //                whichever is larger, (NO_OF_MOLSPECS elts)
	// packetMWDs     are each MWD tree
    // E              is empty space at the end, trying to prevent having to repeat any
    //                communication steps if the mwds exceed the size of the packet

	(*outStatePacket)->stateTooBig = False;
	(*outStatePacket)->noMoreReactions = state.noMoreReactions;
	(*outStatePacket)->time = state.time;
    (*outStatePacket)->deltatemp = state.deltatemp;
	(*outStatePacket)->globalAllMonomer = monomerCount() + getConvertedMonomer();

    // start of array containing total number of each species
	pcount *molCounts = (pcount*)(*outStatePacket+1);

	// copy particle counts
	memcpy(molCounts,state.ms_cnts, sizeof(pcount)*NO_OF_MOLSPECS);

	// copy maximum chain lengths
	int *maxChainLens = (int*)(molCounts + NO_OF_MOLSPECS);
	int bytesForMWDs = getMaxChainLens(maxChainLens);

    int requiredSize = stateCommHeaderBytes() + bytesForMWDs;
//	RANK printf("Required size = %d = %d + %d bytes\n",requiredSize,stateCommHeaderBytes(),bytesForMWDs);
		
    if (requiredSize > stateCommSize) {
		printf("Node %d: state too big for packet (detected)\n", myid);
		printf("Needed %d but have only %d\n", requiredSize, stateCommSize);
		(*outStatePacket)->stateTooBig = True;
		exit(EXIT_FAILURE); // This should not happen
	}
	RANK printf("communication space efficiency = %.0f %%\n", (float)requiredSize/(float)stateCommSize * 100);

	int maxChainLensSize = 0; // size required to store max len of every species in sys
	for (int s=0; s < NO_OF_MOLSPECS; s++) {
		maxChainLensSize += state.arms[s];
	}
	pcount *packetMwds = (pcount*)(maxChainLens + maxChainLensSize);

	// copy mwds into packet
	int pos = 0;
	for (int i=0; i < NO_OF_MOLSPECS; i++) {
		if (isSimpleMol(i)) {
			pos++;
		}
		else {
#ifdef EXPLICIT_SYSTEM_STATE
			if (state.arms[i] == 1) {
				memset(packetMwds,0,maxChainLens[pos]*sizeof(pcount));
				for (int m=0; m<state.ms_cnts[i]; m++) {
					chainLen len = state.expMols[i].mols[m];
					packetMwds[len-1]++;
				}
				packetMwds += maxChainLens[pos];
				pos++;
			}
			else {
				// Hack alert: for data generated for publication, complex species are stored 
				// explicitely in the paket.
				chainLen *compMols = state.expMols[i].mols;
				int bytes = state.ms_cnts[i]*sizeof(chainLen)*state.arms[i];
				memcpy(packetMwds,compMols,bytes);
				maxChainLens[pos] = state.ms_cnts[i]; // for complex species, this stores the total number of molecules, NOT the amx chain length
				IFVERBOSE {
					printf("Before stirr: %lld molecules of %s in packet: ",state.ms_cnts[i],name(i));
					for (int t=0; t<state.ms_cnts[i]*state.arms[i]; t++) {
						printf(" %d", *((chainLen*)packetMwds+t));
					}
					printf("\n");
				}
				packetMwds = (pcount*)((char*)packetMwds + bytes);
				pos++;

				for (int a=1; a<state.arms[i]; a++) {
					maxChainLens[pos++] = 0;
				}
			}
#else
			for (int a=0; a<state.arms[i]; a++) {
				pcount *leaves = state.mwds[i][a].mwd_tree+state.mwds[i][a].maxEntries-1;
				memcpy(packetMwds,leaves,maxChainLens[pos]*sizeof(pcount));
				packetMwds += maxChainLens[pos];
				pos++;
			}
#endif
		}
	}

	if (checkPointing) {
//		char fname[MAX_FILENAME_LEN] = "out.\0";
//		strAppend(fname, ".");
// 
// work in progress
//
		FILE *chk = fopen("out.chk", "w");
		fwrite(*outStatePacket, 1, stateCommSize, chk);
		fclose(chk);
	}

}

int comparitorComplex(const void *m1_, const void *m2_) {
	chainLen *m1 = (chainLen*)m1_;
	chainLen *m2 = (chainLen*)m2_;

	for (int i=0; i<currentComparisonComplexity; i++) {
		if (m1[i] > m2[i])
			return 1;
		if (m1[i] < m2[i])
			return (-1);	
	}
	return 0;
}

/* Takes data from a packet and places it into the current state. The total
 * particles for each species is calculated from the final MWD in the case of
 * polymeric species.
 */
void commToState(StatePacket **inStatePacket) {
	int nextStateCommSizeBase = 0;

	printf("commToState running...\n");

	if ((*inStatePacket)->noMoreReactions == numprocs) {
		// All nodes have no events possible
		printf("No events possible, exiting.\n");
		exit(EXIT_FAILURE);
	}
	
    state.time = (*inStatePacket)->time/(double)numprocs;
	state.deltatemp = (*inStatePacket)->deltatemp/(double)numprocs;

	pcount *speciesCounts = (pcount*)((*inStatePacket) + 1);
	int *maxChainLens = (int*)(speciesCounts + NO_OF_MOLSPECS);	

    // size required to store max len of every species in sys
	int maxChainLensSize = 0;
	for (int s=0; s < NO_OF_MOLSPECS; s++) {
		maxChainLensSize += state.arms[s];
	}

	pcount *mwd = (pcount*)(maxChainLens + maxChainLensSize);
	IFVERBOSE printf("maxChainLensSize = %d\n",maxChainLensSize);
	int pos=0;

	// copy mwd trees into state and take some particles
	for (int i=0; i < NO_OF_MOLSPECS; i++) {
		if (isSimpleMol(i)) {
			// divide particles up amongst processes
			state.ms_cnts[i] = takeSome(speciesCounts[i]);
			pos++;
			continue;
		}

	
		state.ms_cnts[i] = 0;
		
#ifdef EXPLICIT_SYSTEM_STATE
		if (state.arms[i] == 1) { // is poly species

			int maxChainLen = maxChainLens[pos];
			IFVERBOSE printf("maxChainLen[%s] = %d\n",name(i),maxChainLen);

			pcount molsAddedSoFar = 0;
			for (int j=0; j < maxChainLen; j++) {
				pcount mols = takeSome(*(mwd++));
//				printf("mols = %lld\n",mols);
				state.ms_cnts[i] += mols;
				chainLen length = j+1;
				while (mols > 0) {
					if (molsAddedSoFar >= state.expMols[i].maxMolecules) {
						IFVERBOSE printf("masf = %lld mm = %lld\n",molsAddedSoFar,state.expMols[i].maxMolecules);
						printf("system state expanded for %s in commToState\n", name(i));
						state.expMols[i].maxMolecules *= 1.5;
						chainLen *old = state.expMols[i].mols;
						state.expMols[i].mols = malloc(sizeof(chainLen)*state.expMols[i].maxMolecules);
						memcpy(state.expMols[i].mols, old, sizeof(chainLen)*molsAddedSoFar);
						free(old);
					}
					state.expMols[i].mols[molsAddedSoFar++] = length;
					mols--;
				}
			}

			nextStateCommSizeBase += maxChainLen*sizeof(pcount);
			pos++;

		}
		else { // is a complex species
			int start = leaveSome(maxChainLens[pos]);
			pcount molecules = takeSome(maxChainLens[pos]);
			IFVERBOSE printf("After stirr: %lld molecules of %s\n", molecules,name(i));

			int end  = start + molecules;
			IFVERBOSE printf("start = %d end = %d\n",start,end);
			IFVERBOSE printf("maxChainLens[pos] = %d\n",maxChainLens[pos]);
			chainLen *list = (chainLen*)mwd;

			// The complex species must be sorted! This is because different processes
			// will merge them into the packet in different orders resulting in problems
			// when the different processes then take molecules to resume simulation with.
			currentComparisonComplexity = state.arms[i];
			qsort(list, maxChainLens[pos], sizeof(chainLen)*state.arms[i],comparitorComplex);

			state.ms_cnts[i] = 0;
			for (int c=start; c<end; c++) {
				chainLen *mol = list + state.arms[i]*c;
				state.ms_cnts[i]++;
				adjustMolCnt_order1(i, mol, state.arms[i], 1);
			}
			
			IFVERBOSE {
				printf("lengths of %s in packet: ",name(i));
				for (int t=0; t<maxChainLens[pos]*state.arms[i]; t++) {
					printf(" %d", list[t]);
				}
				printf("\n");
			}
			
			mwd = (pcount*)((char*)mwd + sizeof(chainLen)*maxChainLens[pos]*state.arms[i]);

			pos += state.arms[i];
			IFVERBOSE {
				printf("After stirr2: %lld molecules of %s: ", state.ms_cnts[i],name(i));
				for (int t=0; t < state.ms_cnts[i]*state.arms[i]; t++) {
					printf(" %d",state.expMols[i].mols[t]);
				}
				printf("\n");
			}
			
		}
#else
		pcount ms_cnts_tmp_min = state.localParticles; // for some reason LONG_LONG_MAX is not defined
		pcount ms_cnts_tmp_max = 0;
	
		for (int a=0; a<state.arms[i]; a++) {

			while (maxChainLens[pos] > state.mwds[i][a].maxEntries) {
				fprintf(stderr,"commToState: expanding tree memory for species %s(arm=%d)\n",name(i),a);
				free(state.mwds[i][a].mwd_tree);
				state.mwds[i][a].maxEntries *= 2;
				state.mwds[i][a].mwd_tree = malloc(sizeof(pcount)*(2*state.mwds[i][a].maxEntries-1));
			}
			memset(state.mwds[i][a].mwd_tree,'\0',(2*state.mwds[i][a].maxEntries-1)*sizeof(pcount));

			// new MWD becomes leaves in probability tree
			int entries = state.mwds[i][a].maxEntries;
			pcount *leaves = state.mwds[i][a].mwd_tree+entries-1;

			int maxChainLen = maxChainLens[pos];
			assert (maxChainLen <= entries);

			pcount ms_cnts_tmp   = 0;
			for (int j=0; j < maxChainLen; j++) {
				leaves[j] = takeSome(*(mwd++));
				ms_cnts_tmp += leaves[j];
			}
			/* When the total number of arms of complex species are divided
			 * up, some complex species may get "torn apart" because of the 
			 * way the takeSome() function works. For this reason, the 
			 * ms_cnts is set to the lowest number of whole molecules.
			 */
			if (ms_cnts_tmp < ms_cnts_tmp_min)
				ms_cnts_tmp_min = ms_cnts_tmp;
			if (ms_cnts_tmp > ms_cnts_tmp_max)
				ms_cnts_tmp_max = ms_cnts_tmp;

			nextStateCommSizeBase += maxChainLen*sizeof(pcount);
			buildTreeFromLeaves(state.mwds[i][a].mwd_tree,leaves,entries);
			pos++;
		}
		state.ms_cnts[i] = ms_cnts_tmp_min;
		if (ms_cnts_tmp_max - ms_cnts_tmp_min > 0)
			printf("Warning: %lld molecules of species %s temporarily lost upon communication\n", ms_cnts_tmp_max - ms_cnts_tmp_min, name(i));
#endif
	}

	//calculate size of next state comm
	nextStateCommSizeBase += stateCommHeaderBytes();
	stateCommSize = (float)nextStateCommSizeBase*PACKET_SIZE_FACTOR;
	printf("next stateCommSizeBase = %d\n",nextStateCommSizeBase);

	REACTION_PROBABILITY_TREE_INIT
	
#ifdef MONO_AUDIT
	// Global monomer audit
	RANK { 
		if ((*inStatePacket)->globalAllMonomer != GLOBAL_MONOMER_PARTICLES) {
			printf("Warning: global monomer audit has failed. Discrepancy = %lld\n", (*inStatePacket)->globalAllMonomer - GLOBAL_MONOMER_PARTICLES);
		}
		else {
			printf("Global monomer audit passed.\n");
		}
	}
	// setup data for local monomer audits during the next simulation phase
	state.currentMonomerMolecules = monomerCount();
	state.initialMonomerMolecules = state.currentMonomerMolecules + getConvertedMonomer();
	IFVERBOSE printf("imm = %lld\n",state.initialMonomerMolecules);
#endif
}

void checkMWDs(pcount *mwd, int totalLength, char *str) {
	for (int i=0; i<totalLength; i++) {
		if (mwd[i] > state.localMonomerParticles)
			printf("******************* %s: Bad data after stirring totalLen = %d\n",str, totalLength);
	}
}

/* Combines the state information which has been sent. Cyber 
 * equivalent of stirring.
 */
void stirr(void *in_, void *inout_, int *len, MPI_Datatype *datatype) {

    StatePacket *in    = (StatePacket*)in_;
    StatePacket *inout = (StatePacket*)inout_;
// *****************************************************
//	return;
// *****************************************************

	if (in->stateTooBig || inout->stateTooBig) {
		inout->stateTooBig = True;
		printf("Error on node %d: state too big (carried)\n", myid);
		exit(EXIT_FAILURE); // Should not happen if we set the stateCommSize a bit larger than needed
	}

	inout->time += in->time;
	inout->deltatemp += in->deltatemp;
	inout->globalAllMonomer += in->globalAllMonomer;

	pcount *specCounts_in = (pcount*)(in+1);
	int    *maxChainLens_in = (int*)(specCounts_in+NO_OF_MOLSPECS);

	int maxChainLensSize = 0; // size required to store max len of every species in sys
	for (int s=0; s < NO_OF_MOLSPECS; s++) {
		maxChainLensSize += state.arms[s];
	}

	pcount *mwd_in = (pcount*)(maxChainLens_in+maxChainLensSize);
	pcount *specCounts_inout = (pcount*)(inout+1);
	int    *maxChainLens_inout = (int*)(specCounts_inout+NO_OF_MOLSPECS);
	pcount *mwd_inout = (pcount*)(maxChainLens_inout+maxChainLensSize);
	pcount *mwd_inout_start = mwd_inout;

	pcount *mwd_tmp = malloc(stateCommSize);
	pcount *mwd_pos = mwd_tmp;
	int totalLength = 0, miscBytes = 0;

	// merge species counts
	int pos = 0;
	for (int i=0; i<NO_OF_MOLSPECS; i++) {

		specCounts_inout[i] += specCounts_in[i];
		if (isSimpleMol(i)) {
			maxChainLens_inout[i] = 0;
			pos++;
			continue;
		}

#ifdef EXPLICIT_SYSTEM_STATE
		if (i >= MAXPOLY) {
			chainLen *mwd_pos_old = (chainLen*)mwd_pos;

			int inoutBytes = maxChainLens_inout[pos]*sizeof(chainLen)*state.arms[i];
			memcpy(mwd_pos,mwd_inout,inoutBytes);
			mwd_pos = (pcount*)((char*)mwd_pos + inoutBytes);

			int inBytes = maxChainLens_in[pos]*sizeof(chainLen)*state.arms[i];
			memcpy(mwd_pos,mwd_in,inBytes);
			mwd_pos = (pcount*)((char*)mwd_pos + inBytes);

			maxChainLens_inout[pos] += maxChainLens_in[pos];

			IFVERBOSE {
				printf("stirr: %d molecules of %s in stirred packet: ", maxChainLens_inout[pos], name(i));
				for (int t=0; t<maxChainLens_inout[pos]*state.arms[i]; t++) {
					printf(" %d", mwd_pos_old[t]);
				}
				printf("\n");
			}

			miscBytes += inoutBytes + inBytes;
			pos += state.arms[i];
			mwd_inout = (pcount*)((char*)mwd_inout + inoutBytes);
			mwd_in    = (pcount*)((char*)mwd_in + inBytes);
		}
		else {
			/* Allow for the situation which will arise in which a mwd tree in some node
			 * has undergone a size increase before the corresponding mwd in another node.
			 */
			int common = min_value(maxChainLens_in[pos],maxChainLens_inout[pos]);
			checkMWDs(mwd_in,common,"mwd_in");
			checkMWDs(mwd_inout,common,"mwd_inout");
			dvecAdd(mwd_pos,mwd_in,mwd_inout,common);

			assert(state.arms[i] == 1);
			
			if (maxChainLens_in[pos] > maxChainLens_inout[pos]) {
				for (int j=common; j < maxChainLens_in[pos]; j++) {
					mwd_pos[j] = mwd_in[j];
				}
			}
			else {
				for (int j=common; j < maxChainLens_inout[pos]; j++) {
					mwd_pos[j] = mwd_inout[j];
				}
			}
			mwd_in += maxChainLens_in[pos];
			mwd_inout += maxChainLens_inout[pos];
			maxChainLens_inout[pos] = max_value(maxChainLens_inout[pos], maxChainLens_in[pos]);
			checkMWDs(mwd_pos,maxChainLens_inout[pos],"mwd_pos");
			mwd_pos += maxChainLens_inout[pos];
			totalLength += maxChainLens_inout[pos];
			pos++;
		}
#else
		for (int a=0; a<state.arms[i]; a++) {

			/* Allow for the situation which will arise in which a mwd tree in some node
			 * has undergone a size increase before the corresponding mwd in another node.
			 */
			int common = min_value(maxChainLens_in[pos],maxChainLens_inout[pos]);
			checkMWDs(mwd_in,common,"mwd_in");
			checkMWDs(mwd_inout,common,"mwd_inout");
			dvecAdd(mwd_pos,mwd_in,mwd_inout,common);
			
			if (maxChainLens_in[pos] > maxChainLens_inout[pos]) {
				for (int j=common; j < maxChainLens_in[pos]; j++) {
					mwd_pos[j] = mwd_in[j];
				}
			}
			else {
				for (int j=common; j < maxChainLens_inout[pos]; j++) {
					mwd_pos[j] = mwd_inout[j];
				}
			}
			mwd_in += maxChainLens_in[pos];
			mwd_inout += maxChainLens_inout[pos];
			maxChainLens_inout[pos] = max_value(maxChainLens_inout[pos], maxChainLens_in[pos]);
			checkMWDs(mwd_pos,maxChainLens_inout[pos],"mwd_pos");
			mwd_pos += maxChainLens_inout[pos];
			totalLength += maxChainLens_inout[pos];
			pos++;
		}
#endif
	}

	mwd_inout = (pcount*)(maxChainLens_inout+maxChainLensSize);

	int newBytesRequiredForComm = stateCommHeaderBytes() + sizeof(pcount)*totalLength + miscBytes;

	if (newBytesRequiredForComm > stateCommSize) {
		inout->stateTooBig = True;
		printf("warning: state too big (stirr), not yet dealt with\n");
		exit(EXIT_FAILURE);
	}
	else {
		memcpy(mwd_inout_start,mwd_tmp,sizeof(pcount)*totalLength + miscBytes);
		checkMWDs(mwd_tmp,totalLength,"mwd_tmp");
		free(mwd_tmp);
	}
}

// ************* end communication code **********


int MPI_Allreduce_wrapper(void *sendbuf, void *recvbuf, int count,
								  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm ) {
	reduces++;
	Timer t;
    startTimer(&t);
    int r = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);

	int64_t rtime = readTimer(&t);
	total_rtime += rtime;
	
    printf("(T) reduce time = %I64d us for %d bytes\n", rtime, stateCommSize);
	lastReduceTime = rtime;
    return r;
						
}

int compute() {

	RANK file_write_state(START);

/* End if all of the following criteria are met:
    - max simtime is exceeeded
    - max number of events is exceeeded
    - max conversion number is exceeeded 
   In any case, always end when the max walltime has been exceeeded */
	while (((state.time < MAX_SIM_TIME * 60)
		|| (state.events < MAX_EVENTS) 
		|| (state.conversion < MAX_CONVERSION))
		&& (readTimerSec(&state.wallTime) < MAX_WALL_TIME * 60))
	{

		Timer work;
		startTimer(&work);
		unsigned long long prev_events = state.events;

		react();

		int64_t wtime = readTimer(&work);
		total_wtime += wtime;
		printf("total walltime = %I64d\n",total_wtime);
		state.currentMonomerMolecules = monomerCount();
		state.conversion = conversion();

#ifdef MONO_AUDIT
		monomerAudit("post reactions");
#endif

		printf("work time on node %d = %I64d us\n", myid, wtime);
		RANK printf("total number of events on node %d = %lld\n", myid, state.events);
		RANK printf("computational speed on node %d = %.4f events/us\n", myid, (float)(state.events-prev_events)/(float)wtime);
		printf("Pre stirr conversion on node %d = %.8f\n", myid, state.conversion);
		RANK printf("about to synch, time = %f\n", state.time);

#ifndef NO_COMM
		int reduceRes = 0;

		StatePacket *outStatePacket = NULL;
		StatePacket *inStatePacket = NULL;
	
		MPI_Op myOp;
		MPI_Datatype state_t;
		
		/* All packets transmitted must be of the same type and size,
		 * but each process does not know what size the other processes
		 * require to send their system state. So, first an adequate size 
		 * is agreed on.
		 */

		int reqSize_send = requiredStateCommSize();
		int reqSize_received;

		reduceRes = MPI_Allreduce(&reqSize_send, &reqSize_received, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		printf("Suggested size is %d, used is %d\n", reqSize_send, reqSize_received);
		stateCommSize = reqSize_received * 1.2; // Some overprovisioning

		// Next, the packet is created and send

		free (outStatePacket);
		free (inStatePacket);

		outStatePacket = (StatePacket*)malloc(stateCommSize);
		inStatePacket = (StatePacket*)malloc(stateCommSize);

		// create packet to send, and allocate memory for received packet
		stateToComm(&outStatePacket,&inStatePacket);

		// Explain to MPI how type is defined.
		// This must be redone each time since each packet is a single
		// value whose size changes with each send.
		MPI_Type_contiguous(stateCommSize, MPI_CHAR, &state_t); 
		MPI_Type_commit( &state_t ); 
		MPI_Op_create( stirr, True, &myOp ); 
			
		// perform reduction
		reduceRes = MPI_Allreduce_wrapper(outStatePacket, inStatePacket, 1, state_t, myOp, MPI_COMM_WORLD);

		// free operation
		MPI_Op_free(&myOp);

		commToState(&inStatePacket); // update the state after stirring
#ifdef MONO_AUDIT
		monomerAudit("post communicate");
#endif


#endif
	
#ifdef CHANGE_SEED
		// generate new seed
		int nextSeed = wtime+lastReduceTime+fetchpid();
		dsfmt_gv_init_gen_rand(nextSeed);
//		printf("new random seed is %d\n",nextSeed);
#endif

		// Recalculate what's needed
		state.temp = state.basetemp + state.deltatemp;
		state.conversion = conversion();
#ifdef CALCFREEVOLUME
		state.freeVolumeFraction = (VF0 + ALPHA_P * (state.temp - TG_P)) * state.conversion + (VF0 + ALPHA_M * (state.temp - TG_M)) * (1 - state.conversion);
#endif
#ifdef CALCMOMENTSOFDIST
		state.zerothMoment = 1;
		state.firstMoment = 1;
		state.secondMoment= 1;

		for (int i = 0; i < NO_OF_MOLSPECS; i++) {

			if (i >= MAXSIMPLE) {
				state.zerothMoment += state.ms_cnts[i];
				int arms = state.arms[i];

				for (int j = 0; j < state.ms_cnts[i]; j++) {
					chainLen *lengths = &state.expMols[i].mols[arms*j];
					chainLen totalLen = 0;

					for (int a = 0; a < arms; a++) {
						totalLen += lengths[a];
					}
					state.firstMoment += totalLen;
					state.secondMoment += totalLen * totalLen;
				}
			}
		}
#endif

		// Print to synced results to screen and file
		RANK print_state_summary();
		RANK file_write_state(PROFILES);

		state.nextSynchTime += state.synchTime;
		state.nextSynchEvents += (int)((float)state.synchEvents/(float)numprocs);

		printf("Post stirr conversion on node %d = %.8f\n", myid, state.conversion);
    }
   
    return 0;
}

int main(int argc, char *argv[]) {
	RANK printf("Simply version 0.98 beta, released on 2016-06-06\n");
	RANK printf("Program compiled at %s on %s\n",__TIME__,__DATE__);
    system("echo \"starting up on host $HOSTNAME\"");

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	startTimer(&state.wallTime);
	
	printf("numprocs = %d, myid = %d\n", numprocs, myid);

    state.time = 0;// from skeleton, possibly remove in real parallel program?

	// setup defaults that might get overwritten by user options
	state.synchTime = SYNCH_TIME_INTERVAL;
	state.synchEvents = SYNCH_EVENTS_INTERVAL;

#ifndef SEED
    int seed = fetchpid();
#else
	int seed = SEED + myid;
#endif

	// parse command line args:
	char *help = "\n\
-s <rnd seed>\n\
-e <synch events interval>\n\
-t <synch time interval>\n\
-h help\n";
#if defined(_MSC_VER)
	if (argc >= 2) {
		printf("Command line arguments are currently not supported under Windows");
		exit(EXIT_FAILURE);
	}
#elif defined(__GNUC__)
	char *opts = "n:s:t:e:hd";
	for(int c = getopt(argc,argv,opts); c != -1; c = getopt(argc,argv,opts)) {
		switch (c) {
			case -1:
				printf("This should never happen: getopt = -1 not caught\n");
				exit (-1);
			case 's':
				printf ("usr arg: seed = %s\n",optarg);
				sscanf(optarg,"%d",&seed);
				break;
			case 'h':
				printf("Possible arguments:%s",help);
				exit(EXIT_SUCCESS);
				break;
			case 'e':
				printf("synch events = %s\n",optarg);
				sscanf(optarg,"%lld",&state.synchEvents);
				break;
			case 't':
				printf("synch time is every %s seconds\n",optarg);
				sscanf(optarg,"%lf",&state.synchTime);
				break;
		}
	}
#endif

	printf("seed =  %d\n",seed);
    initSysState(seed);

	RANK print_kinetic_model();
//  if (myid == HOST_ID)  while (1) { } //  so I can attach to a debugger when running as multiple processes
  
    printf("Local particles = %lld my seed = %d\n", state.localMonomerParticles, seed);
    RANK print_state();

	// Run the simulation
    compute();

    printf("total time: chatting = %I64d ms avg chat = %I64d working = %I64d\n",total_rtime,(reduces>0?(total_rtime/reduces):0),total_wtime);
	printf("parallel efficiency = %.1lf\n",(float)total_wtime/(float)(total_wtime+total_rtime)*100);

    RANK printf("events = %lld\n",state.events);
    RANK printf("DONE: time = %lf\n",state.time);
	RANK printf("(T) Wall time = %.2f seconds\n", readTimerSec(&state.wallTime));

#ifdef EXPLICIT_SYSTEM_STATE
	RANK compressState();
#endif

	RANK file_write_state(FULL);
	RANK printf("Conversion = %f\n", state.conversion);
    MPI_Finalize();

	exit(EXIT_SUCCESS);
}


