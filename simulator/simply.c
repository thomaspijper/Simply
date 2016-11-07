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

// Block compilers that do not identify as GCC or MSVC
#if defined(__GNUC__)
  // Pass
#elif defined(_MSC_VER)
  #if _MSC_VER < 1900
    #error Visual Studio versions older than VS 2015 are not supported.
  #endif
#else
  #error The compiler you are trying to use is not supported.
#endif

// Change POSIX C SOURCE version for pure C99 and C11 compilers
#if defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
  #if !defined(_POSIX_C_SOURCE) || _POSIX_C_SOURCE < 200112L
    #undef _POSIX_C_SOURCE
    #define _POSIX_C_SOURCE 200112L
  #endif
#endif
 
// Headers usable under both Windows and Linux
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <inttypes.h>
#include <assert.h>
#include <limits.h>
#include <float.h>

// MPI header
#if defined(__GNUC__)
  #include <mpi.h>
#elif defined(_MSC_VER)
  #include "mpi.h"
#endif

// Miscellaneous headers
#if defined(__GNUC__)
  #include <unistd.h>	// POSIX headers, getpid()
  #define _GNU_SOURCE
  #define _BSD_SOURCE
#elif defined(_MSC_VER)
// Compilation succeeds even without these #includes...?
  //#include <process.h>
  //#define WIN32_LEAN_AND_MEAN
  //#include <windows.h>
#endif

// Our own headers
#include "dSFMT.h"
#include "argparse.h"
#include "minunit.h"
#include "genpolymer.h"
#include "simply.h"

// Debugging functions
//#define DEBUG 1
//#define TRACE 1
//#define NO_COMM 1
#define MONO_AUDIT
#define EXPLICIT_SYSTEM_STATE

// Control of random number generation
//#define CHANGE_SEED 1

//#define SCALING 1

#define VERBOSE 1
#define IFVERBOSE if (VERBOSE == 1)
#define IFVERBOSELONG if (VERBOSE == 2)
#define RANK if (myid == 0)
#define VERBOSENODES 1
#define IFVERBOSENODES if (VERBOSENODES == 1)

// Aliases used for printing to files
#define FULL 0
#define PROFILES 1
#define START 2

// Aliases used for printing to screen
#define PRESTIRR 0
#define POSTSTIRR 1

#define HOST_ID 0
#define START_MWD_SIZE 512 // must be a power of 2
#define INIT_STATE_COMM_SIZE (6 * sizeof(pcount) * START_MWD_SIZE)
#define TREE_SAMPLE	10
#define MAX_FILENAME_LEN 100
#define MAX_FILE_SIZE 1048576

#define AVOGADRO 6.022140857E23

/* Number of PRNs to generate at a time.
   Should not be smaller than 382 and must be an even number.
   Changing this value will lead to different PRNs being generated. */
#define PRNG_ARRAY_SIZE 1000

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

unsigned long long total_wtime = 0, total_rtime = 0, reduces = 0;

const int startsize = START_MWD_SIZE;

int myid = -1;

void print_kinetic_model(void);
void monomerAudit(const char *str);
void print_state(void);
       
INLINE pcount max_value(pcount x, pcount y) {
	return (x < y ? y : x);
}

INLINE pcount min_value(pcount x, pcount y) {
	return (x > y ? y : x);
}

// Variables needed for pseudorandom number generation
static w128_t dummy0[PRNG_ARRAY_SIZE / 2 + 1];
static w128_t dummy1[PRNG_ARRAY_SIZE / 2 + 1];
static w128_t dummy2[PRNG_ARRAY_SIZE / 2 + 1];
static double *rndArray0 = (double *)dummy0;
static double *rndArray1 = (double *)dummy1;
static double *rndArray2 = (double *)dummy2;
static size_t rndCounter0 = PRNG_ARRAY_SIZE;
static size_t rndCounter1 = PRNG_ARRAY_SIZE;
static size_t rndCounter2 = PRNG_ARRAY_SIZE;

/* Returns a pseudorandom number.
     randomProb(0) returns number in interval (0,1)
     randomProb(1) returns number in interval [0,1)
     randomProb(2) returns number in interval (0,1]
     randomProb(3) returns number in interval [0,1] -- currently (0,1)
     randomProb(4) forces all arrays to be recomputed
 */
INLINE static double randomProb(int x) {

	if ((x == 0)
		|| (x == 3)) { // interval [0,1] is not available with dSFMT.h, so we use (0,1) instead
		if (rndCounter0 == PRNG_ARRAY_SIZE) {
			dsfmt_gv_fill_array_open_open(rndArray0, PRNG_ARRAY_SIZE); // interval (0,1)
			rndCounter0 = 0;
		}
		double rnd = rndArray0[rndCounter0];
		rndCounter0 += 1;
		return rnd;
	}
	else if (x == 1) {
		if (rndCounter1 == PRNG_ARRAY_SIZE) {
			dsfmt_gv_fill_array_close_open(rndArray1, PRNG_ARRAY_SIZE); // interval [0,1)
			rndCounter1 = 0;
		}
		double rnd = rndArray1[rndCounter1];
		rndCounter1 += 1;
		return rnd;
	}
	else if (x == 2) {
		if (rndCounter2 == PRNG_ARRAY_SIZE) {
			dsfmt_gv_fill_array_open_close(rndArray2, PRNG_ARRAY_SIZE); // interval (0,1]
			rndCounter2 = 0;
		}
		double rnd = rndArray2[rndCounter2];
		rndCounter2 += 1;
		return rnd;
	}
	else if (x == 4) { // Force array to be recalculated
		rndCounter0 = PRNG_ARRAY_SIZE;
		rndCounter1 = PRNG_ARRAY_SIZE;
		rndCounter2 = PRNG_ARRAY_SIZE;
		return 0;
	}
	else {
		printf("You've reached an unavailable option for randomProb. Exiting...\n");
		exit(EXIT_FAILURE);
	}
}

/* Makes each process take the correct number of particles, ensuring conservation of
   particles accross all systems.
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
    tgt_ix = 2 * p - 1;
  }
  *arr = new_arr;
  free (old_arr);
}

void startTimer(Timer *t) {
#if defined(__GNUC__)
	gettimeofday(&t->start_t,NULL);
#elif defined(_MSC_VER)
	QueryPerformanceCounter(&t->start_t);
#endif
}

uint64_t readTimer(Timer *t) {
#if defined(__GNUC__)
	gettimeofday(&t->end_t,NULL);
	return (1000000 * t->end_t.tv_sec + t->end_t.tv_usec) - (1000000 * t->start_t.tv_sec + t->start_t.tv_usec);
#elif defined(_MSC_VER)
	QueryPerformanceCounter(&t->end_t);
	long long elapsedTicks = (t->end_t.QuadPart - t->start_t.QuadPart);
	return ((elapsedTicks * 1000000) / frequency.QuadPart);
#endif
}

double readTimerSec(Timer *t) {
#if defined(__GNUC__)
	gettimeofday(&t->end_t,NULL);
	double begin = (double)t->start_t.tv_sec + t->start_t.tv_usec / 1000000.0;
	double end = (double)t->end_t.tv_sec + t->end_t.tv_usec / 1000000.0;
	return (end - begin);
#elif defined(_MSC_VER)
	QueryPerformanceCounter(&t->end_t);
	QueryPerformanceCounter(&t->end_t);
	long long elapsedTicks = (t->end_t.QuadPart - t->start_t.QuadPart);
	return ((double)elapsedTicks / (double)frequency.QuadPart);
#endif
}

void dumpTree(mwdStore *x) {
	for (int i = 1; i <= x->maxEntries; i *= 2) {
		for (int j = 0; j < i; j++) {
			printf (" %lld",x->mwd_tree[i - 1 + j]);
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
}

/* Add 'diff' molecules of chain length (leave_ind + 1) of type
   spec_ind into tree. Increase tree size if necessary.
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
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		int arms = state.arms[i];
		if (i >= MAXSIMPLE) {
			printf("compressing %lld molecules of %s \n", state.ms_cnts[i], name(i));
			// find maximum chain length to determine space required
			chainLen maxLen = 0;
			for (int j = 0; j < state.ms_cnts[i]; j++) {
				chainLen *lengths = &state.expMols[i].mols[arms*j];
				chainLen totalLen = 0;

				for (int a = 0; a < arms; a++) {
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
			size_t bytes = (2 * state.mwds[i][0].maxEntries - 1) * sizeof(pcount);
			free(state.mwds[i][0].mwd_tree);
			state.mwds[i][0].mwd_tree = malloc(bytes);
			memset(state.mwds[i][0].mwd_tree, 0, bytes);

			// build up MWDs
			for (int j = 0; j < state.ms_cnts[i]; j++) {
				chainLen *lengths = &state.expMols[i].mols[arms*j];
				chainLen totalLen = 0;

				for (int a = 0; a < arms; a++) {
					totalLen += lengths[a];
				}

				IFVERBOSELONG printf("totalLen = %d maxEntries = %d\n", totalLen, state.mwds[i][0].maxEntries);
				int leavesOffset = state.mwds[i][0].maxEntries - 2;
				state.mwds[i][0].mwd_tree[leavesOffset + totalLen]++;
			}
		}
	}
}

/*
 * Converts compact representation to explicit representation. Requires debugging.
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

	return cnt;
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
#ifdef SIMULATEHEATING // Use of TREE_UPDATE_BODY isn't sufficient, so we'll update the entire tree (fast and simple)
	REACTION_PROBABILITY_TREE_INIT
#else
	switch (prevReact)
		TREE_UPDATE_BODY
#endif
}

INLINE double toConc(long long ps) {
	return (ps / AVOGADRO / state.volume);
}

/* Appends string s2 to string s1. Both strings should be defined.
*/
void strAppend(char *s1, const char *s2) {
	size_t len1 = strlen(s1);
	size_t len2 = strlen(s2);
	int count = (int)min_value(len2, MAX_FILENAME_LEN - len1 - 1);

	// Concatenate the strings in a safe way
  #if defined(_MSC_VER)
	int r = strncat_s(s1, MAX_FILENAME_LEN, s2, count);
	if (r != 0) {
		printf("Error while creating file name\nError code reported by strncat_s is %d\n", r);
		exit(EXIT_FAILURE);
	}
  #elif defined(__GNUC__)
	strncat(s1, s2, count);
	s1[MAX_FILENAME_LEN - 1] = "\0";
  #endif
}

/* Opens a file
*/
void fileOpen(FILE **stream, const char *filename, const char *mode) {
  #if defined(_MSC_VER)
	int r = fopen_s(stream, filename, mode);
	if (r != 0) {
		printf("Could not open file %s\nError code produced by fopen_s is %d\n", filename, r);
		exit(EXIT_FAILURE);
	}
  #elif defined(__GNUC__)
	*stream = fopen(filename, mode);
	if (*stream == NULL) {
		printf("Could not open file %s\n", filename);
		exit(EXIT_FAILURE);
	}
  #endif
}

/* Writes state information to files.
*/
void file_write_state(int mode) {
	int i, j, offset, length;
	char timeStr[MAX_FILENAME_LEN];
	snprintf(timeStr, MAX_FILENAME_LEN, "%d", (int)roundf(state.time));

	// Initialize writing of concentrations
	char concfname[MAX_FILENAME_LEN] = "\0";
	strAppend(concfname, "concentrations");
	strAppend(concfname, ".csv");
	FILE *conc;
	fileOpen(&conc, concfname, "a");

	// Initialize writing of rates
	char ratesfname[MAX_FILENAME_LEN] = "\0";
	strAppend(ratesfname, "rates");
	strAppend(ratesfname, ".csv");
	FILE *rates;
	fileOpen(&rates, ratesfname, "a");

	// Initialize writing of rate coefficients
	char ratecoeffsfname[MAX_FILENAME_LEN] = "\0";
	strAppend(ratecoeffsfname, "ratecoeffs");
	strAppend(ratecoeffsfname, ".csv");
	FILE *ratecoeffs;
	fileOpen(&ratecoeffs, ratecoeffsfname, "a");

	if (mode == START) {
		// Write headers for concentrations
		fprintf(conc, "Simulation time (s);Conversion");
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
		fprintf(rates, "Simulation time (s);Conversion");
#ifdef SIMULATEHEATING
		fprintf(rates, ";Temperature (K)");
#endif
		for (i = 0; i < NO_OF_REACTIONS; i++) {
			fprintf(rates, ";%s (mol L^-1 s^-1)", rname(i));
		}
		fprintf(rates, "\n");

		// Write headers for rate coefficients
		fprintf(ratecoeffs, "Simulation time (s);Conversion");
#ifdef SIMULATEHEATING
		fprintf(ratecoeffs, ";Temperature (K)");
#endif
		for (i = 0; i < NO_OF_REACTIONS; i++) {
			if (state.reactions[i].arg_ms2 == NO_MOL) {
				fprintf(ratecoeffs, ";%s (s^-1)", rname(i));
			}
			else {
				fprintf(ratecoeffs, ";%s (L mol^-1 s^-1)", rname(i));
			}
		}
		fprintf(ratecoeffs, "\n");

	}
	else if (mode == PROFILES) {

		// Write concentrations
		fprintf(conc, "%f;%f", state.time, state.conversion);
#ifdef SIMULATEHEATING
		fprintf(conc, ";%.2f", state.temp);
#endif
		for (i = 0; i < NO_OF_MOLSPECS; i++) {
			fprintf(conc, ";%e", 1e6*toConc(state.ms_cnts[i]));
		}
#ifdef CALCMOMENTSOFDIST
		fprintf(conc, ";%e", 1e6*toConc(state.momentDist[0] - 1));
		fprintf(conc, ";%e", 1e6*toConc(state.momentDist[1] - 1));
		fprintf(conc, ";%e", 1e6*toConc(state.momentDist[2] - 1));
		fprintf(conc, ";%e", ((double)state.momentDist[1] / (double)state.momentDist[0]));
		fprintf(conc, ";%e", ((double)state.momentDist[2] / (double)state.momentDist[1]));
		fprintf(conc, ";%e", ((double)state.momentDist[2] * (double)state.momentDist[0] / ((double)state.momentDist[1] * (double)state.momentDist[1])));
#endif
		fprintf(conc, "\n");

		// Write rates
		fprintf(rates, "%f;%f", (state.time), state.conversion);
#ifdef SIMULATEHEATING
		fprintf(rates, ";%.2f", state.temp);
#endif
		for (i = 0; i < NO_OF_REACTIONS; i++) {
			fprintf(rates, ";%e", state.reactProbTree[i + REACT_PROB_TREE_LEAVES - 1] / AVOGADRO / state.volume);
		}
		fprintf(rates, "\n");

		// Write rates coefficients
		fprintf(ratecoeffs, "%f;%f", (state.time), state.conversion);
#ifdef SIMULATEHEATING
		fprintf(ratecoeffs, ";%.2f", state.temp);
#endif
		for (i = 0; i < NO_OF_REACTIONS; i++) {
			if (state.reactions[i].arg_ms2 == NO_MOL) { // Unimolecular reaction
				fprintf(ratecoeffs, ";%e", state.reactions[i].rc);
			}
			else if (state.reactions[i].arg_ms1 == state.reactions[i].arg_ms2) { // Bimolecular reaction with identical reactants
				fprintf(ratecoeffs, ";%e", (state.reactions[i].rc * state.volume * AVOGADRO / SZYMANSKI));
			}
			else { // Bimolecular reaction with nonidentical reactants
					fprintf(ratecoeffs, ";%e", (state.reactions[i].rc * state.volume * AVOGADRO));
			}
		}
		fprintf(ratecoeffs, "\n");

	}
	else if (mode == FULL) {
	// Write MWD of each poly/complex species

		char id[MAX_FILENAME_LEN];
		snprintf(id, MAX_FILENAME_LEN, "%d", myid);

		for (i = 0; i < NO_OF_MOLSPECS; i++) {
			if ((i >= MAXSIMPLE) && (i < MAXPOLY)) {
				char distfname[MAX_FILENAME_LEN] = "\0";
				strAppend(distfname, name(i));
				strAppend(distfname, "-");
				strAppend(distfname, timeStr);
				strAppend(distfname, "-node");
				strAppend(distfname, id);
				strAppend(distfname, ".csv");
				FILE *dist;
				fileOpen(&dist, distfname, "a");

				fprintf(dist, "Chain length;Particle count;Concentration (umol/L)\n");
				offset = state.mwds[i][0].maxEntries - 1;
				length = 1;
				for (j = offset; j < 2 * offset + 1; j++) {
					if (state.mwds[i][0].mwd_tree[j] > 0) {
						fprintf(dist, "%d;%lld;%e\n", length, state.mwds[i][0].mwd_tree[j], 1e6*toConc(state.mwds[i][0].mwd_tree[j]));
					}
					length++;
				}
				fclose(dist);
			}
		}
	}
	fclose(conc);
	fclose(rates);
	fclose(ratecoeffs);
}

/* Picks a random reaction by scanning over the rates.
*
* Since the old scan results are still present, it is
* not necessary to start recomputing the values from
* index 0. It would be sufficient to start with the
* lowest index of all molecules involved in previous reaction.
* This would probably pay off for systems with a higher number
* of reactions.
*/
#if defined(_MSC_VER)
__forceinline int pickRndReaction();
#elif defined(__GNUC__)
inline int pickRndReaction() __attribute__((always_inline));
#endif

INLINE int pickRndReaction() {
    react_prob	rate;
    int			i;
    react_prob	rnd, origRnd;
	static int	prevReact;

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
		if (rate < 0) {
			printf("Fatal error: rate on node %d is %lf. Exiting...\n\n", myid, rate);
			exit(EXIT_FAILURE);
		}
		else {
			return -1;
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
/*void readDataFile(const char *fname, TimesVals *p) {
	FILE *f;
	fileOpen(&f, fname, "r");

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

}*/

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
		//printf("Slave receiving...\n");
		for(int task = MPI_Recv_int_wrap(); task != SETUP_END; task = MPI_Recv_int_wrap()) {
			if (task == SETUP_CONVDATA) {
				int howmany;
			    MPI_Bcast(&howmany, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(state.timeCalcData.ts, howmany, MPI_FLOAT, 0, MPI_COMM_WORLD);
				MPI_Bcast(state.timeCalcData.xs, howmany, MPI_FLOAT, 0, MPI_COMM_WORLD);
				state.timeCalcData.maxIx = howmany;
			}
			else if (task == SETUP_SYSTEMSCALES) {
				int howmany;
			    MPI_Bcast(&howmany, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(state.sysScaleData.ts, howmany, MPI_FLOAT, 0, MPI_COMM_WORLD);
				MPI_Bcast(state.sysScaleData.xs, howmany, MPI_FLOAT, 0, MPI_COMM_WORLD);
				state.sysScaleData.maxIx = howmany;
				//printf("%d system scales received from Master\n",state.sysScaleData.maxIx);
			}
			else {
				printf ("\nError on node %d: bad setup task\n\n", myid);
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
	state.noMoreReactionsLocal = 0;
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
	state.momentDist[0] = 1;
	state.momentDist[1] = 1;
	state.momentDist[2] = 1;
#endif
#ifdef CALCFREEVOLUME
	/* free volume */
	state.freeVolumeFraction = (VF0 + ALPHA_P * (state.temp - TG_P)) * state.conversion + (VF0 + ALPHA_M * (state.temp - TG_M)) * (1 - state.conversion);
#endif
    /* MWDS initialise those ms_cnts not set to 0. rely on calloc for the rest */
    MWD_INITS

	for (i = 0; i < NO_OF_REACTIONS; i++) {
		state.ms_cnts[i] = takeSome(state.ms_cnts[i]);
	}
	
	/* Quick initialization of reaction counts*/
	for (i = 0; i < NO_OF_MOLSPECS; i++) {
		state.react_cnts[i] = 0;
	}

    /* initialise mwds */
    /* enter all mols into mwds */
	for (i = 0; i < NO_OF_MOLSPECS; i++) {
		for (int a = 0; a < MAX_ARMS; a++) {
	        state.mwds[i][a].maxEntries = startsize;
	        state.mwds[i][a].mwd_tree = calloc(2 * startsize - 1, sizeof(pcount));
		}
    }

	state.initialMonomerMolecules = monomerCount();
	state.currentMonomerMolecules = state.initialMonomerMolecules;
	state.volume = (state.localMonomerParticles/AVOGADRO)/MONOMERCONCENTRATION;

	REACTIONS_INIT

	int tmp[NO_OF_MOLSPECS] = ARMS_INIT;
	memcpy(state.arms, tmp, NO_OF_MOLSPECS*sizeof(int));
	memset(state.reactProbTree, 0, (2*REACT_PROB_TREE_LEAVES-1)*sizeof(probability));
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
		pcount lastIx = state.ms_cnts[spec_index]-1;
		pcount rndIx = (pcount)((double)randomProb(3)*(double)lastIx + 0.5); // Fast rounding by adding 0.5, then casting to int
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

INLINE const char *name(int index) {
    switch (index) 
        MOLECULENAMES
    
    return 0;
}


INLINE const char *rname(int index) {
	switch (index) 
		REACTIONNAMES
	
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
		if (i < MAXSIMPLE) {
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

/* React body
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

    RANK printf("Performing reactions until next sycnchronization point is reached...");
    
	do {
       
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

		// No events are possible
		if (reactionIndex == -1) {
			state.noMoreReactionsLocal = 1;
			break;
		}

		state.react_cnts[reactionIndex] += 1;

		react1_ind = state.reactions[reactionIndex].arg_ms1;

#ifdef EXPLICIT_SYSTEM_STATE
		prm1 = pickRndMolecule_order1(react1_ind, react1_lens, &react1_arms);
#else
		prm1 = pickRndMolecule(react1_ind, react1_lens, &react1_arms);
#endif

		if (state.reactions[reactionIndex].arg_ms2 != NO_MOL) {
			react2_ind = state.reactions[reactionIndex].arg_ms2;
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
		if ((prod1_lens[0] >= CHAINLENLIMIT) ||
			(prod2_lens[0] >= CHAINLENLIMIT)) {
			printf("\nError: a particle of species %s on node %d has reached its maximum chain length.\nWriting data and exiting...\n\n", name(prod1_ind), myid);
			printf("\nNode %d exited with errors\n", myid);
			exit(EXIT_FAILURE);
		}

#ifdef CALCMOMENTSOFDIST
		// Adjust moments for disappearance reactant 1
		if (react1_ind >= MAXSIMPLE) {
			state.momentDist[0] -= 1;
			if (react1_ind < MAXPOLY) {
				state.momentDist[1] -= react1_lens[0];
				state.momentDist[2] -= react1_lens[0] * react1_lens[0];
			}
			else {
				int totalLength = 0;
				for (int arm = 0; arm < react1_arms; arm++) {
					totalLength += react1_lens[arm];
				}
				state.momentDist[1] -= totalLength;
				state.momentDist[2] -= totalLength * totalLength;
			}
			// Adjust moments for disappearance reactant 2
		}
		if (react2_ind >= MAXSIMPLE) {
			state.momentDist[0] -= 1;
			if (react2_ind < MAXPOLY) {
				state.momentDist[1] -= react2_lens[0];
				state.momentDist[2] -= react2_lens[0] * react2_lens[0];
			}
			else {
				int totalLength = 0;
				for (int arm = 0; arm < react2_arms; arm++) {
					totalLength += react2_lens[arm];
				}
				state.momentDist[1] -= totalLength;
				state.momentDist[2] -= totalLength * totalLength;
			}
		}
		// Adjust moments for appearance product  1
		if (prod1_ind >= MAXSIMPLE) {
			state.momentDist[0] += 1;
			if (prod1_ind < MAXPOLY) {
				state.momentDist[1] += prod1_lens[0];
				state.momentDist[2] += prod1_lens[0] * prod1_lens[0];
			}
			else {
				int totalLength = 0;
				for (int arm = 0; arm < prod1_arms; arm++) {
					totalLength += prod1_lens[arm];
					state.momentDist[1] += totalLength;
					state.momentDist[2] += totalLength * totalLength;
				}
			}
		// Adjust moments for appearance product 2
		}
		if (prod2_ind >= MAXSIMPLE) {
			state.momentDist[0] += 1;
			if (prod2_ind < MAXPOLY) {
				state.momentDist[1] += prod2_lens[0];
				state.momentDist[2] += prod2_lens[0] * prod2_lens[0];
			}
			else {
				int totalLength = 0;
				for (int arm = 0; arm < prod2_arms; arm++) {
					totalLength += prod2_lens[arm];
					state.momentDist[1] += totalLength;
					state.momentDist[2] += totalLength * totalLength;
				}
			}
		}
#endif

#ifdef RECALCCONVERSION
		// Blindly recalculating the conversion is cheaper than first checking whether the conversion needs to be recalculated
		state.currentMonomerMolecules = monomerCount();
		state.conversion = conversion();
#endif

#ifdef SIMULATEHEATING
		state.deltatemp += (state.reactions[reactionIndex].energy / state.volume);
		state.temp = state.basetemp + state.deltatemp;
#endif

		probability     rndtime;
		probability     rate;
		ptime			deltatime;

		rndtime = randomProb(2);
		rate = state.reactProbTree[0];

		deltatime = (-log(rndtime)) / rate;
		state.time += deltatime;

#ifdef COOLINGRATE
		if (state.deltatemp != 0.0) { // Don't cool when not possible
			state.deltatemp *= exp(-1 * COOLINGRATE * deltatime);
			state.temp = state.basetemp + state.deltatemp;
		}
#endif

#ifdef CALCFREEVOLUME
		state.freeVolumeFraction = (VF0 + ALPHA_P * (state.temp - TG_P)) * state.conversion + (VF0 + ALPHA_M * (state.temp - TG_M)) * (1 - state.conversion);
#endif

#if defined(DEBUG)
		//	monomerAudit(MONOMER_INDEX);
		assert(rate >= 0);
#endif

		state.events++;
	}

	while (((state.time < state.nextSynchTime) || (state.synchTime == 0)) 
		  &&
		  ((state.events < state.nextSynchEvents) || (state.synchEvents == 0)));

    RANK printf("done!\n");

}

// Prints the mwds and molecule counts
 void print_state() {
    int             i, j, offset;
	printf("\n\n");
    for (i = 0; i < NO_OF_MOLSPECS; i++) {
        if (state.ms_cnts[i] > 0) {
            printf("Species %s (%lld): ", name(i), state.ms_cnts[i]);
            if (i < MAXSIMPLE) {
                printf("%lld\n", state.ms_cnts[i]);
            } else {
				pcount cnts[MAX_ARMS];
				offset = state.mwds[i][0].maxEntries;
				for (j = 0; j < offset; j++) {
					pcount total = 0;
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
	printf("\n");
}

void print_state_summary(int m, ptime *simtimes, float *simconversions, double *simtemps, pcount *statecnts, pcount *statemomentdists) {

	size_t nameLens[NO_OF_MOLSPECS];
	size_t maxNumLen = strlen("0.0000e+00");
	size_t nodeIDLen = strlen("Node XYZ ");
	size_t sum = 0;
	int nodesToPrint = 1;
	int rankLen;
	int distNo = 3;

	IFVERBOSENODES nodesToPrint = numprocs;

	sum += nodeIDLen;
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		nameLens[i] = max_value(strlen(name(i)), maxNumLen) + 1;
		sum += nameLens[i];
	}

	printf("\n");
	for (int i = 0; i < sum; i++) {
		printf("-");
	}
	printf("\n");
	
	// Header
	if (m == PRESTIRR) {
		printf("Pre-stirr state information\n\n");
	}
	else if (m == POSTSTIRR) {
		printf("Post-stirr state information\n\n");
	}

	// Wall time and simulation time
	printf("Wall time (s) = %lf\n", readTimerSec(&state.wallTime));
	if (m == PRESTIRR) {
		for (int i = 0; i < nodesToPrint; i++) {
			printf("Simulation time (s) on node %d = %f\n", i, simtimes[i]);
		}
	}
	else if (m == POSTSTIRR) {
		printf("Simulation time (s) = %f\n", simtimes[0]);
	}

	// Conversion
	if (m == PRESTIRR) {
		for (int i = 0; i < nodesToPrint; i++) {
			printf("Conversion on node %d = %f\n", i, simconversions[i]);
		}
	}
	else if (m == POSTSTIRR) {
		printf("Conversion = %f\n", simconversions[0]);
	}
	printf("\n");

#ifdef SIMULATEHEATING
	// Temperature
	if (m == PRESTIRR) {
		for (int i = 0; i < nodesToPrint; i++) {
			printf("Temperature (K) on node %d = %.2f\n", i, (roundf(simtemps[i] * 100) / 100));
		}
	}
	else if (m == POSTSTIRR) {
		printf("Temperature (K) = %.2f\n", (roundf(simtemps[0] * 100) / 100));
	}
	printf("\n");
#endif

	// Concentrations (and optionally moments of distribution, NACL, WACL, and PDI)
	for (int i = 0; i < nodeIDLen; i++)
		printf(" ");
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		printf("%s", name(i));
		for (int j = 0; j < (nameLens[i] - strlen(name(i))); j++) {
			printf(" ");
		}
	}
	printf("\n");

	for (int i = 0; i < numprocs; i++) {
		printf("Node %d", i);
		if (i == 0) {
			rankLen = 1;
		}
		else {
			rankLen = (int)log10(i) + 1;
		}
		for (int j = 0; j < (nodeIDLen - (5 + rankLen)); j++) {
			printf(" ");
		}
		for (int j = 0; j < NO_OF_MOLSPECS; j++) {
			char tmp[MAX_FILENAME_LEN];
			snprintf(tmp, MAX_FILENAME_LEN, "%.4E", 1e6*toConc(statecnts[j + i * NO_OF_MOLSPECS]));
			printf("%s", tmp);
			for (int k = 0; k < (nameLens[j] - strlen(tmp)); k++)
				printf(" ");
		}
		printf("\n");
	}

#ifdef CALCMOMENTSOFDIST
	printf("\n");
	for (int i = 0; i < nodeIDLen; i++)
		printf(" ");
	printf("0th moment 1st moment 2nd moment\n");
	for (int i = 0; i < nodeIDLen; i++)
		printf(" ");
	printf("of distrib of distrib of distrib NACL       WACL       PDI        \n");

	for (int i = 0; i < numprocs; i++) {
		printf("Node %d", i);
		if (i == 0) {
			rankLen = 1;
		}
		else {
			rankLen = (int)log10(i) + 1;
		}
		for (int j = 0; j < (nodeIDLen - (5 + rankLen)); j++) {
			printf(" ");
		}
		printf("%.4E ", 1e6*toConc((statemomentdists[0 + i * distNo] - 1)));
		printf("%.4E ", 1e6*toConc((statemomentdists[1 + i * distNo] - 1)));
		printf("%.4E ", 1e6*toConc((statemomentdists[2 + i * distNo] - 1)));
		printf("%.4E ", (double)statemomentdists[1 + i * distNo] / (double)statemomentdists[0 + i * distNo]);
		printf("%.4E ", (double)statemomentdists[2 + i * distNo] / (double)statemomentdists[1 + i * distNo]);
		printf("%.4E ", (double)statemomentdists[2 + i * distNo] * (double)statemomentdists[0 + i * distNo] / ((double)statemomentdists[1 + i * distNo] * (double)statemomentdists[1 + i * distNo]));
		printf("\n");
	}
#endif

	for (int j = 0; j<sum; j++)
		printf("-");
	printf("\n");
	printf("\n");

}

 void state_summary(int m) {

	int distNo = 3; // zeroth, first, and second

	// Collect simulation times
	ptime *workerStateTime = NULL;
	RANK workerStateTime = malloc(sizeof(ptime) * numprocs);
	MPI_Gather(&state.time, 1, MPI_DOUBLE, workerStateTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Collect conversions
	float *workerStateConversion = NULL;
	RANK workerStateConversion = malloc(sizeof(float) * numprocs);
	MPI_Gather(&state.conversion, 1, MPI_FLOAT, workerStateConversion, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

	// Collect simulation temperatures
	double *workerStateTemp = NULL;
#ifdef SIMULATEHEATING
	RANK workerStateTemp = malloc(sizeof(double) * numprocs);
	MPI_Gather(&state.temp, 1, MPI_DOUBLE, workerStateTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

	// Collect species counts
	pcount *workerStateCnts = NULL;
	RANK workerStateCnts = malloc(sizeof(pcount) * numprocs * NO_OF_MOLSPECS);
	MPI_Gather(&state.ms_cnts, NO_OF_MOLSPECS, MPI_UNSIGNED_LONG_LONG, workerStateCnts, NO_OF_MOLSPECS, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

	// Collect moments of distribution
	pcount *workerStateMomentDist = NULL;
#ifdef CALCMOMENTSOFDIST
	RANK workerStateMomentDist = malloc(sizeof(pcount) * numprocs * distNo);
	MPI_Gather(&state.momentDist, distNo, MPI_UNSIGNED_LONG_LONG, workerStateMomentDist, distNo, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
#endif

	// Print the state, then free up memory
	RANK print_state_summary(m, workerStateTime, workerStateConversion, workerStateTemp, workerStateCnts, workerStateMomentDist);

	free(workerStateTime);
#ifdef SIMULATEHEATING
	free(workerStateTemp);
#endif
	free(workerStateCnts);
#ifdef CALCMOMENTSOFDIST
	free(workerStateMomentDist);
#endif
}


void printMaxChainLens(int mode, unsigned *workerMaxChainLens) {

	size_t maxNameLen = 10; // 9 characters and 1 space
	size_t nodeIDLen = strlen("Node XYZ ");
	size_t sum = 0;
	int nodesToPrint = 1;
	int rankLen = 1;
	int armStrLen = 1;

	IFVERBOSENODES nodesToPrint = numprocs;

	sum += maxNameLen;
	sum += nodeIDLen * nodesToPrint;

	printf("\n");
	for (int i = 0; i < sum; i++) {
		printf("-");
	}
	printf("\n");

	if (mode == PRESTIRR) {
		printf("Pre-stirr maxChainLen\n\n");
	}
	else if (mode == POSTSTIRR) {
		printf("Post-stirr maxChainLen\n\n");
	}

	for (int i = 0; i < maxNameLen; i++) {
		printf(" ");
	}
	for (int i = 0; i < nodesToPrint; i++) {
		printf("Node %d", i);
		if (i == 0) {
			rankLen = 1;
		}
		else {
			rankLen = (int)log10(i) + 1;
		}
		for (int j = 0; j < (nodeIDLen - 5 - rankLen); j++) {
			printf(" ");
		}
	}
	printf("\n");

	int pntrPos = 0;
	for (int i = 0; i < MAXPOLY; i++) {
		if (i >= MAXSIMPLE) {
			//Print name
			printf("%s", name(i));
			size_t max = (maxNameLen - strlen(name(i)));
			for (size_t j = 0; j < max; j++) {
				printf(" ");
			}
			//Print for each node
			for (int j = 0; j < nodesToPrint; j++) {
				// Compare arms lengths for multiarm species
				unsigned armLen = 0;
				for (int k = 0; k < state.arms[i]; k++) {
					armLen = max(armLen, workerMaxChainLens[pntrPos + j * TOTAL_ARMS]);
					pntrPos += 1;
				}
				// Print length
				printf("%u", armLen);
				if (armLen == 0) {
					armStrLen = 1;
				}
				else {
					armStrLen = (int)log10(armLen) + 1;
				}
				for (int k = 0; k < (nodeIDLen - armStrLen); k++) {
					printf(" ");
				}
				pntrPos -= 1; // We've gone 1 too far
			}
			printf("\n");
		}
		pntrPos += 1;
	}
	for (int i = 0; i < sum; i++) {
		printf("-");
	}
	printf("\n");
	printf("\n");
}

/*
  Communicate the maximum chain length of each species, return bytes 
  needed to fit MWDs in state packet. In the current implementation, only
  poly species require max chain lengths to be stored.
  When in verbose mode, collect and print max chain lengths for all species.
  
  The width of 'chainLen' is variable but the MPI command requires a certain data
  type to be specified. As such, use UINT to ensure we can always store the value.
*/
size_t communicateMaxChainLens(int mode, unsigned *maxChainLens) {
#ifdef EXPLICIT_SYSTEM_STATE  // Currently only compatible with EXPLICIT_SYSTEM_STATE on

	// Find global maxChainLen for each species
	size_t bytesNeeded = 0;
	unsigned commMaxChainLen;
	int position = 0;

	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		if (i >= MAXSIMPLE) {
			unsigned commChainLen = maxChainLens[position];
			if (i >= MAXPOLY) { // Complex species, for which maxChainLens actually represents the number of particles
				MPI_Allreduce(&commChainLen, &commMaxChainLen, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
				bytesNeeded += commMaxChainLen * sizeof(chainLen) * state.arms[i];
			}
			else { // Poly species
				MPI_Allreduce(&commChainLen, &commMaxChainLen, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
				bytesNeeded += commMaxChainLen * sizeof(pcount);
			}
		}
		position += state.arms[i];
	}

	IFVERBOSE {
		// Collect maxChainLens from all nodes (poly species only), and print them
		unsigned *workerMaxChainLens = NULL;
		RANK workerMaxChainLens = malloc(sizeof(unsigned) * numprocs * MAXPOLY);
		MPI_Gather(maxChainLens, MAXPOLY, MPI_UNSIGNED, workerMaxChainLens, MAXPOLY, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

		RANK printMaxChainLens(mode, workerMaxChainLens);

		free(workerMaxChainLens);
	}

	IFVERBOSE RANK printf("communicateMaxChainLens: bytesNeeded = %zu\n", bytesNeeded);
	return bytesNeeded;
#endif
}


void print_state_conc(void) {
	int i, j, offset, length;

	printf ("System state:\n\ttime: %f node volume = %E L\n",state.time,state.volume);

	for (i = 0; i < NO_OF_MOLSPECS; i++) {
		if (state.ms_cnts[i] > 0) {
			printf (" %s %lld = (%.0f umol/L):\n", name(i), state.ms_cnts[i], 1e6*toConc(state.ms_cnts[i]));
			if (i < MAXSIMPLE) {
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
	printf("\tname\ttype\n");
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
        char *typeStr;
		if (i < MAXMONOMER)
			typeStr = "monomer";
		else if (i < MAXSIMPLE)
			typeStr = "simple";
		else if (i < MAXPOLY)
			typeStr = "poly";
		else
			typeStr = "complex";

		printf("%d:\t%s\t%s\n",i,name(i),typeStr);
	}
	printf("\nKinetic model\n");
	for (int i=0; i<NO_OF_REACTIONS; i++) {
		reaction *p = state.reactions+i;
		char *p1 = (char*)((p->res_ms1!=NO_MOL)?name(p->res_ms1):" ");
		char *p2 = (char*)((p->res_ms2!=NO_MOL)?name(p->res_ms2):" ");
		char *r1 = (char*)((p->arg_ms1!=NO_MOL)?name(p->arg_ms1):" ");
		char *r2 = (char*)((p->arg_ms2!=NO_MOL)?name(p->arg_ms2):" ");
		double k = state.reactions[i].rc;
		printf("%d:\t%s\t+\t%s\t-->\t%s\t+\t%s\t k = %.5e\n",i,r1,r2,p1,p2,k);
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
		if (i >= MAXSIMPLE) {
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

void monomerAudit(const char *str) {

    pcount discrepancy = state.currentMonomerMolecules
					   + getConvertedMonomer()
					   - state.initialMonomerMolecules;

	pcount *workerDiscrepancies = NULL;
	RANK workerDiscrepancies = malloc(sizeof(pcount) * numprocs);
	MPI_Gather(&discrepancy, 1, MPI_UNSIGNED_LONG_LONG, workerDiscrepancies, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

	RANK {
		int auditFailFlag = 0;
		for (int i = 0; i < numprocs; i++) {
			if (workerDiscrepancies[i] != 0) {
				printf("Local monomer audit (%s) failed on node %d! Discrepancy = %lld\n", str, i, workerDiscrepancies[i]);
				auditFailFlag = 1;
			}
		}

		if (auditFailFlag == 0) {
			printf("Local monomer audit (%s) passed on all nodes\n", str);
		}
	}

	free(workerDiscrepancies);
}


// ************* begin communication code ********

size_t stateCommSize = INIT_STATE_COMM_SIZE;

/* Performs a+b=c for an array
 */
INLINE void dvecAdd(pcount *c, pcount *a, pcount *b, unsigned n) {
	for (int i = 0; i < n; i++) {
			c[i] = a[i] + b[i];
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
 * Populates list of the maximum chain length of all species. Simple 
 * species have chain length = 0. Returns bytes needed to store all
 * MWDs.
 */
size_t getMaxChainLens(unsigned *max_len_arr_ptr) {
	size_t bytesNeeded = 0, pos = 0;
	
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		if (i < MAXSIMPLE) {
			max_len_arr_ptr[pos] = 0;
			pos++;
		} else {
#ifdef EXPLICIT_SYSTEM_STATE
			if (state.arms[i] == 1) { // For single arm species, store max chain length as usual
				unsigned maxLen = 0;
				for (pcount x = 0; x < state.ms_cnts[i]; x++) {
					unsigned tmp = state.expMols[i].mols[x];
					if (tmp > maxLen)
						maxLen = tmp;
				}
				max_len_arr_ptr[pos] = maxLen;
				RANK IFVERBOSE printf("max_len_arr_ptr[%s] on node %d = %u\n", name(i), myid, maxLen);
				bytesNeeded += max_len_arr_ptr[pos];
				pos++;
			}
			else { // For complex species, 1st entry stores number of molecules, rest are zero
				max_len_arr_ptr[pos] = (unsigned)state.ms_cnts[i];
				RANK IFVERBOSE printf("max_len_arr_ptr[%s] on node %d = %u\n", name(i), myid, max_len_arr_ptr[pos]);
				pos++;
				for (int a = 1; a < state.arms[i]; a++) {
					max_len_arr_ptr[pos++] = 0;
				}
			}
#else
			int currMax;
			for (int a = 0; a < state.arms[i]; a++) {
				int offset = state.mwds[i][a].maxEntries - 1;
				currMax = offset - 1;
				for (int j = offset; j < 2*offset+1; j++) {
					if (state.mwds[i][a].mwd_tree[j] > 0) {
						currMax = j;
					}
				}
				max_len_arr_ptr[pos] = currMax - offset + 1;
				bytesNeeded += max_len_arr_ptr[pos];

				pos++;
			}
#endif
		}
    }
	bytesNeeded *= sizeof(pcount); // Calculate space needed to store MWDs for all species
	RANK IFVERBOSE printf("getMaxChainLens on node %d: bytesNeeded = %zu\n", myid, bytesNeeded);
	return bytesNeeded;
}

/*
 * Calculate number of bytes required for the header of the state packet
 */
size_t stateCommHeaderBytes(void) {

	return (sizeof(StatePacket)
		  + sizeof(pcount)*NO_OF_MOLSPECS // Mol counts
          + sizeof(unsigned)*TOTAL_ARMS);  // Maximum chain lengths
}

size_t requiredStateCommSize(void) {
	// Allocate memory, get maxChainLens
	unsigned *maxChainLens = malloc(sizeof(unsigned) * TOTAL_ARMS);
	getMaxChainLens(maxChainLens);
	
	// communicate MaxChainLens between workers, finds largest numbers, calculate space needed for storing MWDs
	size_t bytesForMWDs = communicateMaxChainLens(PRESTIRR, maxChainLens);

	free(maxChainLens);
	return stateCommHeaderBytes() + bytesForMWDs;
}

void stateToComm(StatePacket **outStatePacket, StatePacket **inStatePacket) {

	// Packet format
	// [ Header                                    | Body                                 ]
	// [ StatePacket    | pcount[]  | int[]        | pcount[] | pcount[] | ... | pcount[] ]
	// [ outStatePacket | molCounts | maxChainLens | packetMWDs ...                    E  ]
	// outStatePacket is statepacket struct
	// molCounts      is number of each species (NO_OF_MOLSPECS elts)
	// maxChainLens   is the max. chain length stored for each dist, or MIN_MWD_SIZE, 
    //                whichever is larger, (NO_OF_MOLSPECS elts)
	// packetMWDs     are each MWD tree
    // E              is empty space at the end, trying to prevent having to repeat any
    //                communication steps if the mwds exceed the size of the packet

	(*outStatePacket)->stateTooBig = False;
	(*outStatePacket)->noMoreReactions = state.noMoreReactionsLocal;
	(*outStatePacket)->time = state.time;
    (*outStatePacket)->deltatemp = state.deltatemp;
	(*outStatePacket)->globalAllMonomer = monomerCount() + getConvertedMonomer();

    // Define start of array containing total number of each species, then copy particle counts
	pcount *molCounts = (pcount*)(*outStatePacket + 1);
	memcpy(molCounts, state.ms_cnts, sizeof(pcount) * NO_OF_MOLSPECS);

	// Define start of array containing maximum chain lengths for each species, then copy maximum chain lengths
	unsigned *maxChainLens = (unsigned*)(molCounts + NO_OF_MOLSPECS);
	size_t bytesForMWDs = getMaxChainLens(maxChainLens); // While copying, also determine space needed to store MWDs on this worker

	// Safety: check if stateCommSize is large enough
    size_t requiredSize = stateCommHeaderBytes() + bytesForMWDs;
    if (requiredSize > stateCommSize) { // This should not happen anymore
		printf("Node %d: state too big for packet (prestirr)\n", myid);
		printf("Needed %zu but have only %zu\n", requiredSize, stateCommSize);
		(*outStatePacket)->stateTooBig = True;
		exit(EXIT_FAILURE);
	}
	RANK printf("communication space efficiency = %.0f %%\n", (float)requiredSize/(float)stateCommSize * 100);

	// Define start of array containing MWDs of each species
	pcount *packetMwds = (pcount*)(maxChainLens + TOTAL_ARMS);

	// Copy MWDs into packet
	int pos = 0;
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		if (i < MAXSIMPLE) {
			pos++;
		}
		else {
#ifdef EXPLICIT_SYSTEM_STATE
			// First set count for all chainlengths to 0, then populate the MWD packet for this species
			if (state.arms[i] == 1) {
				memset(packetMwds, 0, maxChainLens[pos] * sizeof(pcount));
				for (int m = 0; m < state.ms_cnts[i]; m++) {
					chainLen len = state.expMols[i].mols[m];
					packetMwds[len-1]++;
				}
				packetMwds += maxChainLens[pos];
				pos++;
			}
			else {
				// Hack alert: complex species are stored explicitely in the packet
				chainLen *compMols = state.expMols[i].mols;
				size_t bytes = state.ms_cnts[i] * sizeof(chainLen) * state.arms[i];
				memcpy(packetMwds, compMols, bytes);
				// maxChainLens[pos] = state.ms_cnts[i]; // redundant line, this is already done in getMaxChainLen()
				IFVERBOSE {
					printf("Before stirr: %lld molecules of %s in packet: ", state.ms_cnts[i], name(i));
					for (int t = 0; t < state.ms_cnts[i]*state.arms[i]; t++) {
						IFVERBOSELONG printf(" %d\n", *((chainLen*)packetMwds+t));
					}
					printf("\n");
				}
				packetMwds = (pcount*)((char*)packetMwds + bytes);
				pos++;

				for (int a = 1; a < state.arms[i]; a++) {
					maxChainLens[pos++] = 0;
				}
			}
#else
			for (int a = 0; a < state.arms[i]; a++) {
				pcount *leaves = state.mwds[i][a].mwd_tree+state.mwds[i][a].maxEntries-1;
				memcpy(packetMwds, leaves, maxChainLens[pos] * sizeof(pcount));
				packetMwds += maxChainLens[pos];
				pos++;
			}
#endif
		}
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

	RANK printf("commToState running...\n");
	
	state.noMoreReactions = (*inStatePacket)->noMoreReactions;
    state.time = (*inStatePacket)->time / numprocs;
	state.deltatemp = (*inStatePacket)->deltatemp / numprocs;

	pcount *speciesCounts = (pcount*)((*inStatePacket) + 1);
	unsigned *maxChainLens = (unsigned*)(speciesCounts + NO_OF_MOLSPECS);

	pcount *mwd = (pcount*)(maxChainLens + TOTAL_ARMS);

	int pos = 0;

	// copy mwd trees into state and take some particles
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		if (i < MAXSIMPLE) {
			// divide particles up amongst processes
			state.ms_cnts[i] = takeSome(speciesCounts[i]);
			pos++;
			continue;
		}

		state.ms_cnts[i] = 0;
		
#ifdef EXPLICIT_SYSTEM_STATE
		if (state.arms[i] == 1) { // poly species

			unsigned maxChainLen = maxChainLens[pos];
			RANK IFVERBOSE printf("maxChainLen[%s] = %u\n",name(i),maxChainLen);

			int takeInterval = numprocs;
			int takeCounter = myid; // myid functions as offset
			pcount molsAddedSoFar = 0;
			
			for (unsigned j = 0; j < maxChainLen; j++) {
				pcount mols = *(mwd++);
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
					if (takeCounter == 0) { // take the molecule
						state.expMols[i].mols[molsAddedSoFar++] = length;
						state.ms_cnts[i] += 1;
						takeCounter += takeInterval;
					}
					takeCounter -= 1;
					mols--;
				}
			}

			pos++;

		}
		else { // complex species
			pcount start = leaveSome(maxChainLens[pos]);
			pcount molecules = takeSome(maxChainLens[pos]);
			IFVERBOSE printf("After stirr: %lld molecules of %s\n", molecules, name(i));

			pcount end = start + molecules;
			IFVERBOSE printf("start = %lld end = %lld\n", start, end);
			IFVERBOSE printf("maxChainLens[pos] = %d\n", maxChainLens[pos]);
			chainLen *list = (chainLen*)mwd;

			// The complex species must be sorted! This is because different processes
			// will merge them into the packet in different orders resulting in problems
			// when the different processes then take molecules to resume simulation with.
			currentComparisonComplexity = state.arms[i];
			qsort(list, maxChainLens[pos], sizeof(chainLen)*(size_t)state.arms[i], comparitorComplex);

			state.ms_cnts[i] = 0;
			for (int c = start; c < end; c++) {
				chainLen *mol = list + state.arms[i]*c;
				state.ms_cnts[i]++;
				adjustMolCnt_order1(i, mol, state.arms[i], 1);
			}
			
			IFVERBOSE {
				printf("lengths of %s in packet: ",name(i));
				for (int t = 0; t < maxChainLens[pos]*state.arms[i]; t++) {
					IFVERBOSELONG printf(" %d\n", list[t]);
				}
				printf("\n");
			}
			
			mwd = (pcount*)((char*)mwd + sizeof(chainLen)*maxChainLens[pos]*state.arms[i]);

			pos += state.arms[i];
			IFVERBOSE {
				printf("After stirr2: %lld molecules of %s: ", state.ms_cnts[i],name(i));
				for (int t=0; t < state.ms_cnts[i]*state.arms[i]; t++) {
					IFVERBOSELONG printf(" %d\n",state.expMols[i].mols[t]);
				}
				printf("\n");
			}
			
		}

#else
		pcount ms_cnts_tmp_min = state.localParticles; // for some reason LONG_LONG_MAX is not defined
		pcount ms_cnts_tmp_max = 0;
	
		for (int a = 0; a < state.arms[i]; a++) {

			while (maxChainLens[pos] > state.mwds[i][a].maxEntries) {
				fprintf(stderr,"commToState: expanding tree memory for species %s(arm=%d)\n", name(i), a);
				free(state.mwds[i][a].mwd_tree);
				state.mwds[i][a].maxEntries *= 2;
				state.mwds[i][a].mwd_tree = malloc(sizeof(pcount) * (2 * state.mwds[i][a].maxEntries - 1));
			}
			memset(state.mwds[i][a].mwd_tree,'\0',(2*state.mwds[i][a].maxEntries-1)*sizeof(pcount));

			// new MWD becomes leaves in probability tree
			unsigned entries = state.mwds[i][a].maxEntries;
			pcount *leaves = state.mwds[i][a].mwd_tree + entries - 1;

			unsigned maxChainLen = maxChainLens[pos];
			assert (maxChainLen <= entries);

			pcount ms_cnts_tmp = 0;
			for (unsigned j = 0; j < maxChainLen; j++) {
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

			buildTreeFromLeaves(state.mwds[i][a].mwd_tree,leaves,entries);
			pos++;
		}
		state.ms_cnts[i] = ms_cnts_tmp_min;
		if (ms_cnts_tmp_max - ms_cnts_tmp_min > 0)
			printf("Warning: %lld molecules of species %s temporarily lost upon communication\n", ms_cnts_tmp_max - ms_cnts_tmp_min, name(i));
#endif
	}

	size_t bytesForMWDs = getMaxChainLens(maxChainLens);
	IFVERBOSE communicateMaxChainLens(POSTSTIRR, maxChainLens);

	REACTION_PROBABILITY_TREE_INIT
	
#ifdef MONO_AUDIT
	// Global monomer audit
	RANK { 
		if ((*inStatePacket)->globalAllMonomer != GLOBAL_MONOMER_PARTICLES) {
			printf("Global monomer audit has failed! Discrepancy = %lld\n", (*inStatePacket)->globalAllMonomer - GLOBAL_MONOMER_PARTICLES);
		}
		else {
			printf("Global monomer audit passed.\n");
		}
	}
	// setup data for local monomer audits during the next simulation phase
	state.currentMonomerMolecules = monomerCount();
	state.initialMonomerMolecules = state.currentMonomerMolecules + getConvertedMonomer();
	IFVERBOSELONG printf("imm = %lld\n",state.initialMonomerMolecules);
#endif
}

void checkMWDs(pcount *mwd, unsigned totalLength, char *str) {
	for (unsigned i = 0; i < totalLength; i++) {
		if (mwd[i] > state.localMonomerParticles)
			printf("******************* %s: Bad data after stirring totalLen = %d\n", str, totalLength);
	}
}

/* Combines the state information which has been sent. Cyber 
 * equivalent of stirring.
 */
void stirr(void *in_, void *inout_, int *len, MPI_Datatype *datatype) {

    StatePacket *in    = (StatePacket*)in_;
    StatePacket *inout = (StatePacket*)inout_;

	if (in->stateTooBig || inout->stateTooBig) { // This should not happen anymore
		inout->stateTooBig = True;
		printf("Error on node %d: state too big (carried)\n", myid);
		exit(EXIT_FAILURE);
	}

	inout->noMoreReactions += in->noMoreReactions;
	inout->time += in->time;
	inout->deltatemp += in->deltatemp;
	inout->globalAllMonomer += in->globalAllMonomer;

	// Define layout of packet 'in'
	pcount *specCounts_in = (pcount*)(in+1);
	unsigned *maxChainLens_in = (unsigned*)(specCounts_in + NO_OF_MOLSPECS);

	pcount *mwd_in = (pcount*)(maxChainLens_in + TOTAL_ARMS);

	// Define layout of packet 'inout'
	pcount *specCounts_inout = (pcount*)(inout+1);
	unsigned *maxChainLens_inout = (unsigned*)(specCounts_inout + NO_OF_MOLSPECS);
	pcount *mwd_inout = (pcount*)(maxChainLens_inout + TOTAL_ARMS);

	pcount *mwd_inout_start = mwd_inout;
	pcount *mwd_tmp = malloc(stateCommSize);
	pcount *mwd_pos = mwd_tmp;
	size_t totalLength = 0, miscBytes = 0;

	// Merge species counts and MWDs
	int pos = 0;
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {

		specCounts_inout[i] += specCounts_in[i];
		if (i < MAXSIMPLE) {
			maxChainLens_inout[i] = 0;
			pos++;
			continue;
		}

#ifdef EXPLICIT_SYSTEM_STATE
		if (i >= MAXPOLY) { // Complex species
			chainLen *mwd_pos_old = (chainLen*)mwd_pos;

			size_t inoutBytes = maxChainLens_inout[pos] * sizeof(chainLen) * state.arms[i];
			memcpy(mwd_pos, mwd_inout, inoutBytes);
			mwd_pos = (pcount*)((char*)mwd_pos + inoutBytes);

			size_t inBytes = maxChainLens_in[pos] * sizeof(chainLen) * state.arms[i];
			memcpy(mwd_pos, mwd_in, inBytes);
			mwd_pos = (pcount*)((char*)mwd_pos + inBytes);

			maxChainLens_inout[pos] += maxChainLens_in[pos];

			IFVERBOSE {
				printf("stirr: %u molecules of %s in stirred packet: ", maxChainLens_inout[pos], name(i));
				for (unsigned t = 0; t < maxChainLens_inout[pos] * state.arms[i]; t++) {
					IFVERBOSELONG printf(" %d", mwd_pos_old[t]);
				}
				printf("\n");
			}

			miscBytes += inoutBytes + inBytes;
			pos += state.arms[i];
			mwd_inout = (pcount*)((char*)mwd_inout + inoutBytes);
			mwd_in    = (pcount*)((char*)mwd_in + inBytes);
		}
		else { // Poly species
			/* Allow for the situation which will arise in which a mwd tree in some node
			   has undergone a size increase before the corresponding mwd in another node.

               Note: this may no longer be possible now that nodes first agree on an adequate MWD tree size
               Redundant code?
             */
			unsigned common = min_value(maxChainLens_in[pos],maxChainLens_inout[pos]);
			checkMWDs(mwd_in, common, "mwd_in");
			checkMWDs(mwd_inout, common, "mwd_inout");
			dvecAdd(mwd_pos, mwd_in, mwd_inout, common);

			assert(state.arms[i] == 1);
			
			if (maxChainLens_in[pos] > maxChainLens_inout[pos]) {
				for (unsigned j = common; j < maxChainLens_in[pos]; j++) {
					mwd_pos[j] = mwd_in[j];
				}
			}
			else {
				for (unsigned j = common; j < maxChainLens_inout[pos]; j++) {
					mwd_pos[j] = mwd_inout[j];
				}
			}
			mwd_in += maxChainLens_in[pos];
			mwd_inout += maxChainLens_inout[pos];
			maxChainLens_inout[pos] = max_value(maxChainLens_inout[pos], maxChainLens_in[pos]);
			checkMWDs(mwd_pos,maxChainLens_inout[pos],"mwd_pos");
			mwd_pos += maxChainLens_inout[pos];
			totalLength += maxChainLens_inout[pos];
			IFVERBOSE RANK printf("stirr: maxChainLens_inout[pos] for species %d is %u\n", i, maxChainLens_inout[pos]);
			pos++;
		}
#else
		for (int a = 0; a < state.arms[i]; a++) {

			/* Allow for the situation which will arise in which a mwd tree in some node
			 * has undergone a size increase before the corresponding mwd in another node.
			 */
			unsigned common = min_value(maxChainLens_in[pos], maxChainLens_inout[pos]);
			checkMWDs(mwd_in, common, "mwd_in");
			checkMWDs(mwd_inout, common, "mwd_inout");
			dvecAdd(mwd_pos, mwd_in, mwd_inout, common);
			
			if (maxChainLens_in[pos] > maxChainLens_inout[pos]) {
				for (unsigned j = common; j < maxChainLens_in[pos]; j++) {
					mwd_pos[j] = mwd_in[j];
				}
			}
			else {
				for (unsigned j = common; j < maxChainLens_inout[pos]; j++) {
					mwd_pos[j] = mwd_inout[j];
				}
			}
			mwd_in += maxChainLens_in[pos];
			mwd_inout += maxChainLens_inout[pos];
			maxChainLens_inout[pos] = max_value(maxChainLens_inout[pos], maxChainLens_in[pos]);
			checkMWDs(mwd_pos, maxChainLens_inout[pos], "mwd_pos");
			mwd_pos += maxChainLens_inout[pos];
			totalLength += maxChainLens_inout[pos];
			pos++;
		}
#endif
	}

	mwd_inout = (pcount*)(maxChainLens_inout + TOTAL_ARMS);

	size_t newBytesRequiredForComm = stateCommHeaderBytes() + sizeof(pcount) * totalLength + miscBytes;

	if (newBytesRequiredForComm > stateCommSize) { // This should not happen anymore
		inout->stateTooBig = True;
		printf("Error: state too big (while stirring)\n");
		exit(EXIT_FAILURE);
	}
	else {
		memcpy(mwd_inout_start, mwd_tmp, sizeof(pcount) * totalLength + miscBytes);
		checkMWDs(mwd_tmp, totalLength, "mwd_tmp");
		free(mwd_tmp);
	}
}

int MPI_Allreduce_wrapper(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm ) {
	reduces++;
	Timer t;
    startTimer(&t);
    int r = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);

	uint64_t rtime = readTimer(&t);
	total_rtime += rtime;
	
    RANK printf("Reduce time (us) = %I64d for %zu bytes\n", rtime, stateCommSize);
	lastReduceTime = rtime;
    return r;
}

/* Dump state content -- currently not functional yet (data likely needs to be serialized) */
void dumpStatePacket(StatePacket **inStatePacket, int stateCommSize) {

	// Open the file
	char fname[MAX_FILENAME_LEN] = "\0";
	strAppend(fname, "stateDump");
	FILE *dump;
	fileOpen(&dump, fname, "w");

	// Check if file could be opened
	if (dump == NULL) {
		printf("Could not write to file %s\n", fname);
		exit(EXIT_FAILURE);
	}

	// Write to the file
	int r1 = fwrite(*inStatePacket, stateCommSize, 1, dump);
	printf("wrote %d elements out of %d requested\n", r1, stateCommSize);

	//Close the file
	fclose(dump);
}

// ************* end communication code **********


// ************* start testing code *********

/* We're going to test a few functions and constants that allow for testing. Extensive unit testing is not possible
   as the code is not written modular. We'll have to test only what we can. */  

void test_setup() {
	
}

void test_teardown() {

}

/*  Test whether myid equals 0 */
MU_TEST(myid_test) {
	mu_assert_int_eq(0, myid);
}

/*  Test whether whether PRNG_ARRAY_SIZE is set correctly and the PRNG gives the expected numbers (the latter ensures reproducibility of results) */
MU_TEST(dSFMT_test) {
	char value[MAX_FILENAME_LEN];
	snprintf(value, MAX_FILENAME_LEN, "%d", DSFMT_N64);
	char errormsg[MAX_FILENAME_LEN] = "PRNG_ARRAY_SIZE may not be smaller than \0";
	strAppend(errormsg, value);
	mu_assert(PRNG_ARRAY_SIZE >= DSFMT_N64, errormsg);
	mu_assert((PRNG_ARRAY_SIZE % 2) == 0, "PRNG_ARRAY_SIZE must be an even number");

	// 1) Seed the PRNG, 2) test the value, 3) empty the PRNG arrays
	// Test path 0
	dsfmt_gv_init_gen_rand(1);
	mu_assert_double_eq(0.1193544251137, randomProb(0));
	randomProb(4);

	// Test path 1
	dsfmt_gv_init_gen_rand(1);
	mu_assert_double_eq(0.1193544251137, randomProb(1));
	randomProb(4);

	// Test path 2
	dsfmt_gv_init_gen_rand(1);
	mu_assert_double_eq(0.8806455748863, randomProb(2));
	randomProb(4);

	// Test path 3
	dsfmt_gv_init_gen_rand(1);
	mu_assert_double_eq(0.1193544251137, randomProb(3));
	randomProb(4);

	// Reseed the PRNG for future use
	dsfmt_gv_init_gen_rand(1);
}

/*  Test whether START_MWD_SIZE is a power of two */
MU_TEST(START_MWD_SIZE_test) {
	int value = START_MWD_SIZE;
	mu_assert(value > 0, "START_MWD_SIZE is not a power of two");
	int value_test = value & (~value + 1);
	mu_assert(value == value_test, "START_MWD_SIZE is not a power of two");
}

/*  Test the min_value() and max_value() functions */
MU_TEST(minmax_test) {
	mu_assert_int_eq(9223372036854775808, min_value(9223372036854775808, 9223372036854775809));
	mu_assert_int_eq(9223372036854775809, max_value(9223372036854775808, 9223372036854775809));
}

/*  Test the takeSome() function by taking the sum of all takeSome(i) values for various values of i */
MU_TEST(takeSome_test) {
	for (int i = 0; i < (2 * numprocs); i++) {
		pcount sendValue = takeSome(i);
		pcount receiveValue;
		MPI_Allreduce(&sendValue, &receiveValue, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
		RANK mu_assert_int_eq(i, receiveValue);
	}
}

/*  Test the leaveSome() function by taking the sum of all takeSome(i) values for various values of i, then checking if the pattern in the sums is correct */
MU_TEST(leaveSome_test) {
	pcount *sumvalues = NULL;
	RANK sumvalues = malloc(sizeof(int) * numprocs * 2);
	for (int i = 0; i < (2 * numprocs); i++) { // fetch sum of all leaveSome() values
		pcount sendValue = leaveSome(i);
		pcount receiveValue;
		MPI_Allreduce(&sendValue, &receiveValue, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
		RANK sumvalues[i] = receiveValue;
	}

	int expectedValue = numprocs - 1;
	RANK for (int i = 1; i < (2 * numprocs); i++) {
		if (expectedValue < 0) {
			expectedValue = numprocs - 1;
		}
		mu_assert_int_eq(expectedValue, (sumvalues[i] - sumvalues[i - 1]));
		expectedValue--;
	}
	RANK free(sumvalues);
}

/*  Test the double_arraysize() function by expanding the size of a dummy array and verifying its contents */
MU_TEST(double_arraysize_test) {
	// Create test array
	int curr_size = 16;
	pcount testnumber = 1;
	pcount *arr = malloc(sizeof(pcount) * curr_size);
	for (int i = 0; i < curr_size; i++) {
		arr[i] = i + testnumber;
	}

	// Double array size
	double_arraysize(&arr, curr_size);

	// Evaluate resulting array
	mu_assert_int_eq(1, arr[0]); // First value in array
	int base = 0;
	int period = 2;
	int expectedValue = testnumber;
	int i = 1;
	while (i < curr_size) {
		for (int j = 0; j < period; j++) {
			if (i <= (base + period / 2)) {
				mu_assert_int_eq(expectedValue, arr[i]);
				expectedValue++;
			} else {
				mu_assert_int_eq(0, arr[i]);
			}
			i++;
		}
		base += period;
		period *= 2;
	}

	mu_assert_int_eq(expectedValue, arr[i]); // Last value in array
}

/*  Test the double_arraysize() function by summing the elements of two dummy arrays */
MU_TEST(dvecAdd_test) {
	// Initialize dummy arrays.
	pcount array1[3] = { 0, 1, 2 };
	pcount array2[3] = { 9223372036854775808, 9223372036854775808, 9223372036854775808 };
	pcount array12[3];

	// Perform function
	dvecAdd(array12, array1, array2, 3);

	// Test whether dvecAdd() worked as intended (it should be able to handle ULL sized values)
	for (int i = 0; i < 3; i++) {
		mu_assert_int_eq(9223372036854775808 + i, array12[i]);
	}
}

/* Tests to run */
MU_TEST_SUITE(test_suite) {

	MU_SUITE_CONFIGURE(&test_setup, &test_teardown);

	RANK MU_RUN_TEST(myid_test);
	RANK MU_RUN_TEST(dSFMT_test);
	RANK MU_RUN_TEST(START_MWD_SIZE_test);
	RANK MU_RUN_TEST(minmax_test);
	RANK MU_RUN_TEST(double_arraysize_test);
	RANK MU_RUN_TEST(dvecAdd_test);

	if (numprocs > 1) {
		MU_RUN_TEST(takeSome_test);
		MU_RUN_TEST(leaveSome_test);
	}
}

// ************* end testing code *********


/*  Parse and apply command line options */
void argumentParsing(int argc, char *argv[], int *seed) {

	int arg_seed = 0;
	int arg_syncevents = 0;
	int arg_synctime = 0;
	struct argparse_option options[] = {
		OPT_HELP(),
		OPT_GROUP("Basic options"),
		OPT_INTEGER('s', "seed", &arg_seed, "seed for PRNG"),
		OPT_INTEGER('e', "syncevents", &arg_syncevents, "number of events between each two synchronizations"),
		OPT_INTEGER('t', "time", &arg_synctime, "simulation time between each two synchronizations"),
		OPT_END(),
	};
	struct argparse argparse;
	argparse_init(&argparse, options, usages, 0);
	argparse_describe(&argparse, "\nSimply - kinetic Monte Carlo simulator for polymerizations.", "\n");
	argc = argparse_parse(&argparse, argc, argv);
	if (arg_seed != 0) {
		*seed = arg_seed + myid;
		RANK printf("Overriding compiled seed with user specified seed: %d\n", *seed);
	}
	else {
		RANK printf("Seed: %d\n", *seed);
	}
	if (arg_syncevents != 0) {
		state.synchEvents = arg_syncevents;
		RANK printf("Overriding compiled number of events synchronization interval with user specified number: %lld\n", state.synchEvents);
	}
	else {
		RANK printf("Number of events between each two synchronizations: %lld\n", state.synchEvents);
	}
	if (arg_synctime != 0) {
		state.synchTime = arg_synctime;
		RANK printf("Overriding compiled simulation time synchronization interval with user specified time: %f\n", state.synchTime);
	}
	else {
		RANK printf("Simulation time between each two synchronizations: %f\n", state.synchTime);
	}
}

// ************* end command line parsing *********

int compute(void) {

	RANK file_write_state(START);
	RANK file_write_state(PROFILES);

#ifndef NO_COMM
	StatePacket *outStatePacket = NULL;
	StatePacket *inStatePacket = NULL;
#endif

/* Stop if all of the following criteria are met:
    - max simtime is exceeeded
    - max number of events is exceeeded
    - max conversion number is exceeeded 
   In addition, always stop when the maximum
   walltime has been exceeeded or no more reactions
   are possible on one or more nodes */
	while (((state.time < MAX_SIM_TIME)
		|| (state.events < MAX_EVENTS) 
		|| (state.conversion < MAX_CONVERSION))
		&& (readTimerSec(&state.wallTime) < MAX_WALL_TIME * 60)
		&& (state.noMoreReactions == 0))
	{

		Timer work;
		startTimer(&work);
		unsigned long long prev_events = state.events;

		react();

		int64_t wtime = readTimer(&work);
		total_wtime += wtime;
		RANK printf("Total calculation walltime (us) = %I64d\n", total_wtime);
		state.currentMonomerMolecules = monomerCount();
		state.conversion = conversion();

#ifdef MONO_AUDIT
		monomerAudit("post reactions");
#endif

		RANK printf("Work time (us) on node %d = %I64d\n", myid, wtime);
		RANK printf("Total number of events on node %d = %lld\n", myid, state.events);
		RANK printf("Computational speed on node %d = %.4f events/us\n", myid, (float)(state.events-prev_events)/(float)wtime);
		
		state_summary(PRESTIRR);

		RANK printf("About to synch, time = %f\n", state.time);

#ifndef NO_COMM
		int reduceRes = 0;

		MPI_Op myOp;
		MPI_Datatype state_t;
		
		/* All packets transmitted must be of the same type and size,
		 * but each process does not know what size the other processes
		 * require to send their system state. So, first an adequate size 
		 * is agreed on.
		 */
		stateCommSize = requiredStateCommSize();

		// Discard previous sync information
		if (outStatePacket) {
			free(outStatePacket);
		}
		if (inStatePacket) {
			free(inStatePacket);
		}

        // Create packet to send, and allocate memory for received packet
		outStatePacket = (StatePacket*)malloc(stateCommSize);
		inStatePacket = (StatePacket*)malloc(stateCommSize);

		stateToComm(&outStatePacket,&inStatePacket);

		// Explain to MPI how type is defined.
		// This must be redone each time since each packet is a single
		// value whose size changes with each send.
		MPI_Type_contiguous(stateCommSize, MPI_CHAR, &state_t); 
		MPI_Type_commit(&state_t); 
		MPI_Op_create(stirr, True, &myOp);
			
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
		state.momentDist[0] = 1;
		state.momentDist[1] = 1;
		state.momentDist[2]= 1;

		for (int i = 0; i < NO_OF_MOLSPECS; i++) {

			if (i >= MAXSIMPLE) {
				state.momentDist[0] += state.ms_cnts[i];
				int arms = state.arms[i];

				for (int j = 0; j < state.ms_cnts[i]; j++) {
					chainLen *lengths = &state.expMols[i].mols[arms*j];
					chainLen totalLen = 0;

					for (int a = 0; a < arms; a++) {
						totalLen += lengths[a];
					}
					state.momentDist[1] += totalLen;
					state.momentDist[2] += totalLen * totalLen;
				}
			}
		}
#endif

		// Print synced results to screen and file
		state_summary(POSTSTIRR);
		RANK file_write_state(PROFILES);

		state.nextSynchTime += state.synchTime;
		state.nextSynchEvents += (int)((float)state.synchEvents/(float)numprocs);
    }
   
    return 0;
}

int main(int argc, char *argv[]) {
    system("echo \"Starting up on host $HOSTNAME\"");

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#if defined(_MSC_VER)
	QueryPerformanceFrequency(&frequency);
#endif
	startTimer(&state.wallTime);

	RANK printf("\n");
	RANK printf("Simply version 0.99 beta prerelease\n");
	RANK printf("Program compiled at %s on %s\n",__TIME__,__DATE__);	
	RANK printf("\n");
	RANK printf("Number of nodes = %d\n", numprocs);
	RANK printf("\n");

	// Unit testing
	MU_RUN_SUITE(test_suite);
	MU_REPORT();
	
	state.synchTime = SYNCH_TIME_INTERVAL;
	state.synchEvents = SYNCH_EVENTS_INTERVAL;
#ifndef SEED
    int seed = fetchpid();
#else
	int seed = SEED + myid;
#endif

	argumentParsing(argc, argv, &seed);

	// Initialize the calculation
	RANK printf("\n");
    initSysState(seed);
	RANK print_kinetic_model();
    RANK print_state();

	// Run the simulation
    compute();

	RANK if (state.noMoreReactions != 0) {
		printf("Exiting prematurely as %d node(s) ran out of reactions to perform\n\n", state.noMoreReactions);
	}

    RANK printf("Total time (us): chatting = %I64d (avg = %I64d), working = %I64d\n", total_rtime, (reduces>0?(total_rtime/reduces):0), total_wtime);
	RANK printf("Parallel efficiency = %.1lf\n", (float)total_wtime/(float)(total_wtime+total_rtime)*100);
	
	RANK printf("\n");
    RANK printf("Events on node %d = %lld\n", myid, state.events);
	RANK for (int i = 0; i < NO_OF_REACTIONS; i++) {
		printf("Events for reaction %s on node %d = %lld\n", rname(i), myid, state.react_cnts[i]);
	}
	RANK printf("\n");

	RANK printf("Final conversion = %f\n", state.conversion);
	RANK printf("Final simulation time = %lf\n", state.time);
	RANK printf("Wall time (s) = %.2f\n", readTimerSec(&state.wallTime));

#ifdef EXPLICIT_SYSTEM_STATE
	RANK compressState();
#endif

	RANK file_write_state(FULL);

	MPI_Finalize();

	exit(EXIT_SUCCESS);
}


