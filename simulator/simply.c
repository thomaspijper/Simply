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
  #if __GNUC__ < 5 || (__GNUC__ == 4 && __GNUC_MINOR__ < 3)
    #error GCC versions older than 4.3.x are not supported.
  #endif
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
#include <assert.h>
#include <limits.h>
#include <malloc.h>

// MPI header
#if defined(__GNUC__)
  #include <mpi.h>
#elif defined(_MSC_VER)
  #include "mpi.h"
#endif

// Miscellaneous headers
#if defined(__GNUC__)
  #include <unistd.h>	// POSIX headers, getpid(), gethostname()
  #include <sys/time.h> // gettimeofday()
  #include <errno.h>    // errno
  #include <dirent.h>   // DIR
  #define _GNU_SOURCE
  #define _BSD_SOURCE
#elif defined(_MSC_VER)
  #include <process.h>  // _getpid()
  #define getpid _getpid
#endif

// Our own headers
#include "dSFMT.h"
#include "argparse.h"
#include "minunit.h"
#include "genpolymer.h"
#include "simply.h"

/* The number of PRNs to generate at a time.
Should not be smaller than 382 and must be an even number.
Changing this value will lead to different PRNs being generated. */
#define PRNG_ARRAY_SIZE 10000

// Verbose options
#define VERBOSELEVEL 0
#define NODEVERBOSELEVEL 1

// Debugging options
#define DEBUGLEVEL 0

// Miscellaneous defines
#define START_MWD_SIZE 512 // must be a power of 2
#define INIT_STATE_COMM_SIZE (6 * sizeof(pcount) * START_MWD_SIZE)
#define MAX_FILENAME_LEN 255
#define MAX_HOSTNAME_LEN 64
#define MAX_FILE_SIZE 1048576
#define AVOGADRO 6.022140857E23
#define RANK if (myid == 0)

// Initialize various variables
enum bools { False = 0, True = 1 };
enum { START, PROFILES };             // for writing state info to files
enum { PRESTIRR, POSTSTIRR };         // for printing state info to the screen
int myid = -1;
int numprocs = -1;
int currentComparisonComplexity = -1;
int mwdsNodesMerged = 1;
unsigned long long total_wtime = 0, total_rtime = 0, reduces = 0, lastReduceTime = 0;
static sysState state;
size_t stateCommSize;
static char dirname[MAX_FILENAME_LEN] = "\0";
static char simname[MAX_FILENAME_LEN] = "\0";
static char hostname[MAX_HOSTNAME_LEN] = "\0";
static double *rndArray0;
static double *rndArray2;
static double *logRndArray2;

// Initialize options
static const _Bool monomeraudit = MONO_AUDIT;
static const _Bool changeseed = CHANGESEED;
static const _Bool explicit = EXPLICITSYSTEM;
static const _Bool calcmoments = CALCMOMENTSOFDIST;
static const _Bool simulateheating = SIMULATEHEATING;
static const _Bool recalcconversion = RECALCCONVERSION;
static const _Bool calcfreevolume = CALCFREEVOLUME;
static const double coolingrate = COOLINGRATE;

// Forced inlines
FORCEINLINE_PRE static void updateTree(int prevReact) FORCEINLINE_POST;
FORCEINLINE_PRE static int pickRndReaction(void) FORCEINLINE_POST;
//FORCEINLINE_PRE static int pickRndChainLen(pcount *mwd_tree, int size) FORCEINLINE_POST;

INLINE static pcount max_value(pcount x, pcount y) {
	return (x < y ? y : x);
}

INLINE static pcount min_value(pcount x, pcount y) {
	return (x > y ? y : x);
}

/* Returns a pseudorandom number.
     randomProb(0) returns number in interval (0,1)
     randomProb(1) returns number in interval [0,1) (disabled)
     randomProb(2) returns number in interval (0,1]
     randomProb(3) returns number in interval [0,1] -- currently (0,1)
     randomProb(4) forces all arrays to be recomputed
 */
INLINE static double randomProb(int x) {
	static size_t rndCounter0 = PRNG_ARRAY_SIZE;
	static size_t rndCounter2 = PRNG_ARRAY_SIZE;
	if ((x == 0) || 
		(x == 3)) { // interval [0,1] is not available with dSFMT.h, so we use (0,1) instead
		if (rndCounter0 == PRNG_ARRAY_SIZE) {
			dsfmt_gv_fill_array_open_open(rndArray0, PRNG_ARRAY_SIZE); // interval (0,1)
			rndCounter0 = 0;
		}
		double rnd = rndArray0[rndCounter0];
		rndCounter0 += 1;
		return rnd;
	}
	/*else if (x == 1) {
		if (rndCounter1 == PRNG_ARRAY_SIZE) {
			dsfmt_gv_fill_array_close_open(rndArray1, PRNG_ARRAY_SIZE); // interval [0,1)
			rndCounter1 = 0;
		}
		double rnd = rndArray1[rndCounter1];
		rndCounter1 += 1;
		return rnd;
	}*/
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
		//rndCounter1 = PRNG_ARRAY_SIZE;
		rndCounter2 = PRNG_ARRAY_SIZE;
		return 0;
	}
	else {
		printf("You've reached an unavailable option for randomProb. Exiting...\n");
		exit(EXIT_FAILURE);
	}
}

/* Produces an array of PRNs on the interval (0,1), then takes the natural logarithm of each.
*  Used for time evolution. Allows for loop vectorization (Visual Studio).
*/
INLINE static double logRndArray(void) {
	static size_t logRndCounter2 = PRNG_ARRAY_SIZE;
	if (logRndCounter2 == PRNG_ARRAY_SIZE) {
		dsfmt_gv_fill_array_open_open(logRndArray2, PRNG_ARRAY_SIZE); // interval (0,1)
		logRndCounter2 = 0;
		for (int i = 0; i < PRNG_ARRAY_SIZE; i++) {
			logRndArray2[i] = log(logRndArray2[i]);
		}
	}
	double logRnd = logRndArray2[logRndCounter2];
	logRndCounter2 += 1;
	return logRnd;
}

/* Makes each process take the correct number of particles, ensuring conservation of
   particles accross all systems.
 */
INLINE static pcount takeSome(pcount total) {
	pcount ans = total/numprocs;
	pcount leftovers = total % (pcount)numprocs;
	if (myid < leftovers) {
		ans++;
	}
	return ans;
}

INLINE static pcount leaveSome(pcount total) {
	pcount rough = total/numprocs;
	pcount leftovers = total % (pcount)numprocs;
	pcount correction = (myid >= leftovers) ? leftovers : myid;
	pcount ans = myid * rough + correction;
	return ans;
}

/* Deals with overflow of array containing mwd tree
 */
static void double_arraysize(pcount **arr, size_t curr_size) {
	size_t newsize = 2 * curr_size;
	pcount *new_arr = calloc(sizeof(pcount) * (2 * newsize - 1), 1);
	pcount *old_arr = *arr;

	size_t p      = 1;  /* size of current level */
	size_t src_ix = 0;  /* src index             */
	size_t tgt_ix = 1;  /* target index          */

	new_arr[0] = old_arr[0];
	while (p < newsize) {
		memcpy(&(new_arr[tgt_ix]), &(old_arr[src_ix]), p * sizeof(pcount));
		src_ix += p;
		p = 2*p;
		tgt_ix = 2 * p - 1;
	}
	*arr = new_arr;
	free (old_arr);
}

/* Appends string s2 to string s1. Both strings should be defined.
*/
static void strAppend(char *s1, const char *s2) {
	size_t len1 = strlen(s1);
	size_t len2 = strlen(s2);
	if ((len1 + len2) > MAX_FILENAME_LEN) {
		printf("\nError: concatenation of the strings \n\n    \"%s\"\n\nand\n\n    \"%s\"\n\nwould result in a string length that exceeds the maximum permitted length of %d characters. Aborting...\n", s1, s2, MAX_FILENAME_LEN);
		exit(EXIT_FAILURE);
	}
	size_t count = (size_t)min_value(len2, MAX_FILENAME_LEN - len1 - 1);

	// Concatenate the strings in a safe way
#if defined(_MSC_VER)
	int r = strncat_s(s1, MAX_FILENAME_LEN, s2, count);
	if (r != 0) {
		printf("Error while creating file name\nError code reported by strncat_s is %d\n", r);
		exit(EXIT_FAILURE);
	}
#elif defined(__GNUC__)
	strncat(s1, s2, count);
	s1[MAX_FILENAME_LEN - 1] = '\0'; // just to be on the safe side
#endif
}

// ************* begin timer code ********

/* Start a timer */
static void startTimer(Timer *t) {
#if defined(__GNUC__)
	gettimeofday(&t->start_t,NULL);
#elif defined(_MSC_VER)
	QueryPerformanceCounter(&t->start_t);
#endif
}

/* Read a timer, return microseconds */
static unsigned long long readTimer(Timer *t) {
#if defined(__GNUC__)
	gettimeofday(&t->end_t,NULL);
	return (1000000 * t->end_t.tv_sec + t->end_t.tv_usec) - (1000000 * t->start_t.tv_sec + t->start_t.tv_usec);
#elif defined(_MSC_VER)
	QueryPerformanceCounter(&t->end_t);
	long long elapsedTicks = (t->end_t.QuadPart - t->start_t.QuadPart);
	return ((elapsedTicks * 1000000) / frequency.QuadPart);
#endif
}

/* Read a timer, return seconds */
static double readTimerSec(Timer *t) {
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

/* Get the current time (currently used for debugging only) */
INLINE static void getSystemTimeString(char *timeStampString) {
#if defined(__GNUC__)
	struct timeval fileTime;
	struct tm systemTime;

	gettimeofday(&fileTime, NULL);
	localtime_r(&fileTime.tv_sec, &systemTime);
	long ms = fileTime.tv_usec / 1000;
	snprintf(timeStampString, 27, "[%04d-%02d-%02d %02d:%02d:%02d:%03ld] ", systemTime.tm_year, systemTime.tm_mon, systemTime.tm_mday,
																		    systemTime.tm_hour, systemTime.tm_min, systemTime.tm_sec,
																		    ms);
#elif defined(_MSC_VER)
	FILETIME fileTime;
	SYSTEMTIME systemTime;

	GetSystemTimeAsFileTime(&fileTime);
	FileTimeToSystemTime(&fileTime, &systemTime);
	snprintf(timeStampString, 27, "[%04d-%02d-%02d %02d:%02d:%02d:%03d] ", systemTime.wYear, systemTime.wMonth, systemTime.wDay,
		                                                                   systemTime.wHour, systemTime.wMinute, systemTime.wSecond,
		                                                                   systemTime.wMilliseconds);
#endif
}

// ************* end timer code ********


/* Gets the host name of the system */
static void retreiveHostname(void) {
#if defined(_MSC_VER)
	unsigned len = MAX_HOSTNAME_LEN - 1;
	BOOL r = GetComputerName(hostname, &len);
	if (r == False) {
		DWORD error = GetLastError();
		printf("\nFunction %s failed (error value = %u), aborting...", __FUNCTION__, error);
		exit(EXIT_FAILURE);
	}
#elif defined(__GNUC__)
	int r = gethostname(hostname, MAX_HOSTNAME_LEN - 1);
	if (r == -1) {
		int error = errno;
		printf("\nFunction %s failed (error value = %d), aborting...", __FUNCTION__, error);
	}
#endif
}

/* Parses a string containing a path */
static void parseDirname(char* path) {
#if defined(_MSC_VER)
	// Check and repair a path that was placed between double quotes and than ended with a '\'
	// such as "C:\Data files\" as the final backslash will be seen as an escape character
	size_t len = strlen(path);
	if (path[len - 1] == '"') {
		path[len - 1] = '\0';
	}
	// Append backslash if missing
	len = strlen(path);
	if (path[len - 1] != '\\') {
		strAppend(path, "\\");
	}
	
	// Check if correct
	DWORD r = GetFileAttributes(path);
	if (r == INVALID_FILE_ATTRIBUTES) {
		RANK printf("\nError: invalid path specified (%s)\n", path);
		exit(EXIT_FAILURE);
	}
	else if (r & FILE_ATTRIBUTE_DIRECTORY) {
		if (r & FILE_ATTRIBUTE_READONLY) {
			RANK printf("\nError: path specified is read-only (%s)\n", path);
			exit(EXIT_FAILURE);
		}
	}
	else {
		RANK printf("\nError: invalid path specified (%s)\n", path);
		exit(EXIT_FAILURE);
	}
#elif defined(__GNUC__)
	// Append slash if missing
	size_t len = strlen(path);
	if (path[len - 1] != '/') {
		strAppend(path, "/");
	}

	// Check if correct
	DIR* dir = opendir(path);
	if (dir) {
		closedir(dir);
		if (access(path, W_OK) != 0) {
			RANK printf("\nError: path specified seems to be read-only (%s)\n", path);
			exit(EXIT_FAILURE);
		}
	}
	else if (ENOENT == errno) {
		RANK printf("\nError: invalid path specified (%s)\n", path);
		exit(EXIT_FAILURE);
	}
	else {
		RANK printf("\nUnknown error when trying to verify path (%s)\n", path);
		exit(EXIT_FAILURE);
	}
#endif
}

static void dumpTree(const mwdStore *x) {
	for (size_t i = 1; i <= x->maxEntries; i *= 2) {
		for (size_t j = 0; j < i; j++) {
			printf(" %llu",x->mwd_tree[i - 1 + j]);
		}
		printf("\n");
	}
}

static void dumpAllTrees(void) {
	for (int i =0; i < NO_OF_MOLSPECS; i++) {
		printf("%s:\n",name(i));
		for (int a=0; a < state.arms[i]; a++) {
			dumpTree(&state.mwds[i][a]);
		}
	}
}

INLINE static void adjustTree(const int spec_ind, mwdStore *mwd, int leave_ind, int diff) {
	while (leave_ind >= mwd->maxEntries) {
		printf("tree data structure overflow: %s leaveIx = %d maxEntries = %d\n", name(spec_ind),leave_ind, mwd->maxEntries);
		double_arraysize(&(mwd->mwd_tree), (size_t)mwd->maxEntries);
		mwd->maxEntries *= 2;
	}

	int ind;
	int levelSize = 1;
	while (levelSize <= mwd->maxEntries) {
		ind = (levelSize - 1) + (levelSize * leave_ind) / mwd->maxEntries;
		mwd->mwd_tree[ind] += diff;
		levelSize *= 2;
	}
}

/* Add 'diff' molecules of chain length (leave_ind + 1) of type
   spec_ind into tree. Increase tree size if necessary.
 */
INLINE static void adjustMolCnt(const int spec_ind, chainLen *lengths, int arms, int diff) {
	if (spec_ind < MAXPOLY) {
		adjustTree(spec_ind, &state.mwds[spec_ind][0], lengths[0]-1, diff);
		if (calcmoments) {
			// Adjust moments of distribution for appearance product
			state.momentDist[0]++;
			state.momentDist[1] += lengths[0];
			state.momentDist[2] += lengths[0] * lengths[0];
		}
	}
	else {
		chainLen comb_length = 0;
		for (int a = 0; a < arms; a++) {
			adjustTree(spec_ind, &state.mwds[spec_ind][a], lengths[a] - 1, diff);
			if (calcmoments) {
				comb_length += lengths[a];
			}
		}
		if (calcmoments) {
			// Adjust moments of distribution for appearance product
			state.momentDist[0]++;
			state.momentDist[1] += comb_length;
			state.momentDist[2] += comb_length * comb_length;
		}
	}
}

/*
 * Converts explicit state representation to compact representation for use 
 * when dumping detailed system contents.
 * Only generates leaves of MWD trees.
 */
static void compressState(void) {
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		int arms = state.arms[i];
		if (i >= MAXSIMPLE) {
			// find maximum chain length to determine space required
			unsigned maxLenLocal = 0; // cannot use the variable 'chainLen' type as the type needs to be known by the MPI operation later
			for (int j = 0; j < state.ms_cnts[i]; j++) {
				chainLen *lengths = &state.expMols[i].mols[arms*j];
				chainLen totalLen = 0;

				for (int a = 0; a < arms; a++) {
					totalLen += lengths[a];
				}
				if (totalLen > maxLenLocal)
					maxLenLocal = totalLen;
			}

			// make all nodes consider the same maximum chain length (necessary if we want to merge MWDs between nodes at a later time)
			unsigned maxLen = 0;
			MPI_Allreduce(&maxLenLocal, &maxLen, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
			
			// increase allocated memory to store system, if required
			while (maxLen > (unsigned)state.mwds[i][0].maxEntries) {
				//printf("%s: mwd tree for %s doubled\n", __FUNCTION__, name(i));
				state.mwds[i][0].maxEntries *= 2;
			}
			size_t bytes = (2 * state.mwds[i][0].maxEntries - 1) * sizeof(pcount);
			free(state.mwds[i][0].mwd_tree);
			state.mwds[i][0].mwd_tree = malloc(bytes);
			memset(state.mwds[i][0].mwd_tree, 0, bytes);

			// build up MWDs
			for (int j = 0; j < state.ms_cnts[i]; j++) {
				chainLen *lengths = &state.expMols[i].mols[arms * j];
				chainLen totalLen = 0;

				for (int a = 0; a < arms; a++) {
					totalLen += lengths[a];
				}

				int leavesOffset = state.mwds[i][0].maxEntries - 2;
				state.mwds[i][0].mwd_tree[leavesOffset + totalLen]++;
			}
		}
	}
}

/*
 * Converts compact representation to explicit representation. Currently nonfunctional; requires review/debugging.
 */
static void decompressState(void) {
	printf("decompressStateRunning... ");
	for (int i=0; i < NO_OF_MOLSPECS; i++) {
		memset(state.expMols[i].mols, 0, sizeof(chainLen) * state.expMols[i].maxMolecules * state.arms[i]);
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
			for (chainLen len = 1; len <= offset; len++) {
				while (state.mwds[i][0].mwd_tree[offset - 1 + len] > 0) {
					state.expMols[i].mols[pos++] = len;
				}
			}
		}
	}
	printf("finished!\n");
}

/* Merges MWDs in compact representation to node 0. For use when dumping system contents */
static void mergeMWDs(void) {
	if (numprocs < 2) { // Nothing to do
		return;
	}

	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		if (i >= MAXSIMPLE) {

			int offset = state.mwds[i][0].maxEntries - 1;
			int treeSize = 2 * state.mwds[i][0].maxEntries - 1;

			// Retreive MWDs from all nodes
			pcount *mwd_tree_combined = NULL;
			RANK mwd_tree_combined = malloc((treeSize) * sizeof(pcount) * numprocs);
			MPI_Gather((&state.mwds[i][0])->mwd_tree, treeSize, MPI_UNSIGNED_LONG_LONG, mwd_tree_combined, treeSize, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

			// Merge MWDs
			RANK for (int j = offset; j < 2 * offset + 1; j++) {
				for (int k = 1; k < numprocs; k++) {
					state.mwds[i][0].mwd_tree[j] += mwd_tree_combined[j + treeSize * k];
				}
			}
			free(mwd_tree_combined);
		}
	}
	mwdsNodesMerged = numprocs;
}

/* For O(1) algorithm, uses large amounts of memory */
INLINE static void adjustMolCnt_order1(const int spec_ind, chainLen *lengths, int arms) {
	pcount molecules = state.ms_cnts[spec_ind];

	if (molecules > state.expMols[spec_ind].maxMolecules) {
		state.expMols[spec_ind].maxMolecules = (pcount)((double)state.expMols[spec_ind].maxMolecules * 1.5);
		RANK printf("Explicit sys state expanded for species %s from %llu to %llu\n", name(spec_ind), (molecules-1), state.expMols[spec_ind].maxMolecules);
		chainLen *old = state.expMols[spec_ind].mols;
		state.expMols[spec_ind].mols = malloc(state.expMols[spec_ind].maxMolecules*arms*sizeof(chainLen));
		memcpy(state.expMols[spec_ind].mols, old, sizeof(chainLen)*(molecules-1)*arms);
		free(old);
	}
	chainLen *mempos = state.expMols[spec_ind].mols + (molecules-1)*arms;
	chainLen comb_length = 0;
	for (int i = 0; i < arms; i++) {
		mempos[i] = lengths[i];
		if (calcmoments) {
			comb_length += mempos[i];
		}
	}
	if (calcmoments) {
		// Adjust moments of distribution for appearance product
		state.momentDist[0]++;
		state.momentDist[1] += comb_length;
		state.momentDist[2] += comb_length * comb_length;
	}
}

INLINE static pcount monomerCount(void) {
	pcount cnt = 0;
	for (int i = 0; i < MAXMONOMER; i++) {
		cnt += state.ms_cnts[i];
	}
	return cnt;
}

INLINE static float conversion(void) {
	float conversion = ((float)state.localMonomerParticles * (float)state.scaleFactor - (float)state.currentMonomerMolecules) / ((float)state.localMonomerParticles*state.scaleFactor);
	return conversion;
}

/*
 *  Returns number of molecules of type spec_ind and length 
 *  (leave_ind +1)
 */
static pcount getMolCnt(int spec_ind, int arm, int leave_ind) {
	int ix = leave_ind+state.mwds[spec_ind][arm].maxEntries-1;
	pcount cnt = state.mwds[spec_ind][arm].mwd_tree[ix];

	return cnt;
}

static void dumpReactProbTree(void) {
	int nextLev = 2;
	for (int j=0; j<2*REACT_PROB_TREE_LEAVES-1; j++) {
		printf (" %.40f",state.reactProbTree[j]);
		if (j == nextLev-2) {
			printf("\n");
			nextLev *= 2;
		}
	}
	printf("\n\n");
}

INLINE static void updateTree(int prevReact) {
	RATES_UPDATE_BODY
	if (simulateheating) { // Use of TREE_UPDATE_BODY currently isn't sufficient, so we'll update the entire tree (fast and simple)
		REACTION_PROBABILITY_TREE_INIT
	}
	else {
		switch (prevReact)
			TREE_UPDATE_BODY
	}
}

INLINE static double toConc(long long ps) {
	return (ps / AVOGADRO / state.volume);
}

/* Opens a file
*/
static void fileOpen(FILE **stream, const char *filename, const char *mode) {
#if defined(_MSC_VER)
	int r = fopen_s(stream, filename, mode);
	if (r != 0) {
		RANK printf("Could not open file %s\nError code produced by fopen_s is %d\n", filename, r);
		exit(EXIT_FAILURE);
	}
#elif defined(__GNUC__)
	*stream = fopen(filename, mode);
	if (*stream == NULL) {
		RANK printf("Could not open file %s\n", filename);
		exit(EXIT_FAILURE);
	}
#endif
}

/* Writes state information to files
*/
static void file_write_state(int mode) {

	// Initialize writing of concentrations
	char concfname[MAX_FILENAME_LEN] = "\0";
	strAppend(concfname, dirname);
	strAppend(concfname, simname);
	strAppend(concfname, "concentrations");
	strAppend(concfname, ".csv");
	FILE *conc;
	fileOpen(&conc, concfname, "a");

	// Initialize writing of rates
	char ratesfname[MAX_FILENAME_LEN] = "\0";
	strAppend(ratesfname, dirname);
	strAppend(ratesfname, simname);
	strAppend(ratesfname, "rates");
	strAppend(ratesfname, ".csv");
	FILE *rates;
	fileOpen(&rates, ratesfname, "a");

	// Initialize writing of rate coefficients
	char ratecoeffsfname[MAX_FILENAME_LEN] = "\0";
	strAppend(ratecoeffsfname, dirname);
	strAppend(ratecoeffsfname, simname);
	strAppend(ratecoeffsfname, "ratecoeffs");
	strAppend(ratecoeffsfname, ".csv");
	FILE *ratecoeffs;
	fileOpen(&ratecoeffs, ratecoeffsfname, "a");

	// Initialize writing of reaction events
	char eventsfname[MAX_FILENAME_LEN] = "\0";
	strAppend(eventsfname, dirname);
	strAppend(eventsfname, simname);
	strAppend(eventsfname, "reactionevents");
	strAppend(eventsfname, ".csv");
	FILE *events;
	fileOpen(&events, eventsfname, "a");

	if (mode == START) {
		// Write headers for concentrations
		fprintf(conc, "Simulation time (s);Conversion");
		if (simulateheating) {
			fprintf(conc, ";Temperature (K)");
		}
		for (int i = 0; i < NO_OF_MOLSPECS; i++) {
			fprintf(conc, ";%s (umol/L)", name(i));
		}
		if (calcmoments) {
			fprintf(conc, ";Zeroth moment of dist (umol/L);First moment of dist (umol/L);Second moment of dist (umol/L);Number-average chain length;Weight-average chain length;Polydispersity index");
		}
		fprintf(conc, "\n");
		
		// Write headers for rates
		fprintf(rates, "Simulation time (s);Conversion");
		if (simulateheating) {
			fprintf(rates, ";Temperature (K)");
		}
		for (int i = 0; i < NO_OF_REACTIONS; i++) {
			fprintf(rates, ";%s (mol L^-1 s^-1)", rname(i));
		}
		fprintf(rates, "\n");

		// Write headers for rate coefficients
		fprintf(ratecoeffs, "Simulation time (s);Conversion");
		if (simulateheating) {
			fprintf(ratecoeffs, ";Temperature (K)");
		}
		for (int i = 0; i < NO_OF_REACTIONS; i++) {
			if (state.reactions[i].arg_ms2 == NO_MOL) {
				fprintf(ratecoeffs, ";%s (s^-1)", rname(i));
			}
			else {
				fprintf(ratecoeffs, ";%s (L mol^-1 s^-1)", rname(i));
			}
		}
		fprintf(ratecoeffs, "\n");

		// Write headers for reaction events
		fprintf(events, "Simulation time (s);Conversion");
		if (simulateheating) {
			fprintf(events, ";Temperature (K)");
		}
		for (int i = 0; i < NO_OF_REACTIONS; i++) {
			fprintf(events, ";%s", rname(i));
		}
		fprintf(events, "\n");
	}

	else if (mode == PROFILES) {
		// Write concentrations
		fprintf(conc, "%f;%f", state.time, state.conversion);
		if (simulateheating) {
			fprintf(conc, ";%.2f", state.temp);
		}
		for (int i = 0; i < NO_OF_MOLSPECS; i++) {
			fprintf(conc, ";%e", 1e6*toConc(state.ms_cnts[i]));
		}
		if (calcmoments) {
			fprintf(conc, ";%e", 1e6*toConc(state.momentDist[0] - 1));
			fprintf(conc, ";%e", 1e6*toConc(state.momentDist[1] - 1));
			fprintf(conc, ";%e", 1e6*toConc(state.momentDist[2] - 1));
			fprintf(conc, ";%e", ((double)state.momentDist[1] / (double)state.momentDist[0]));
			fprintf(conc, ";%e", ((double)state.momentDist[2] / (double)state.momentDist[1]));
			fprintf(conc, ";%e", ((double)state.momentDist[2] * (double)state.momentDist[0] / ((double)state.momentDist[1] * (double)state.momentDist[1])));
		}
		fprintf(conc, "\n");

		// Write rates
		fprintf(rates, "%f;%f", state.time, state.conversion);
		if (simulateheating) {
			fprintf(rates, ";%.2f", state.temp);
		}
		for (int i = 0; i < NO_OF_REACTIONS; i++) {
			fprintf(rates, ";%e", state.reactProbTree[i + REACT_PROB_TREE_LEAVES - 1] / AVOGADRO / state.volume);
		}
		fprintf(rates, "\n");

		// Write rates coefficients
		fprintf(ratecoeffs, "%f;%f", state.time, state.conversion);
		if (simulateheating) {
			fprintf(ratecoeffs, ";%.2f", state.temp);
		}
		for (int i = 0; i < NO_OF_REACTIONS; i++) {
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

		// Write number of reaction events
		fprintf(events, "%f;%f", state.time, state.conversion);
		if (simulateheating) {
			fprintf(events, ";%.2f", state.temp);
		}
		for (int i = 0; i < NO_OF_REACTIONS; i++) {
			fprintf(events, ";%llu", state.react_cnts[i]);
		}
		fprintf(events, "\n");

	}

	fclose(conc);
	fclose(rates);
	fclose(ratecoeffs);
	fclose(events);
}

/* Write MWD of each poly/complex species to files
*/
static void file_write_MWDs(void) {
	char timeStr[MAX_FILENAME_LEN];
	snprintf(timeStr, MAX_FILENAME_LEN, "%d", (int)round(state.time));
	char id[MAX_FILENAME_LEN];
	snprintf(id, MAX_FILENAME_LEN, "%d", myid);

	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		if ((i >= MAXSIMPLE) && (i < MAXPOLY)) {
			char distfname[MAX_FILENAME_LEN] = "\0";
			strAppend(distfname, dirname);
			strAppend(distfname, simname);
			strAppend(distfname, name(i));
			strAppend(distfname, "-");
			strAppend(distfname, timeStr);
			strAppend(distfname, ".csv");
			FILE *dist;
			fileOpen(&dist, distfname, "a");

			fprintf(dist, "Chain length;Particle count;Concentration (umol/L)\n");
			int offset = state.mwds[i][0].maxEntries - 1;
			int length = 1;
			for (int j = offset; j < 2 * offset + 1; j++) {
				if (state.mwds[i][0].mwd_tree[j] > 0) {
					fprintf(dist, "%d;%llu;%e\n", length, state.mwds[i][0].mwd_tree[j], (1e6 * toConc(state.mwds[i][0].mwd_tree[j]) / mwdsNodesMerged));
				}
				length++;
			}
			fclose(dist);
		}
	}
}

/* Writes debug information to a dedicated file
*/
static void file_write_debug(const char *str) {

	if (DEBUGLEVEL < 1) { // Do nothing
		return;
	}

	char id[MAX_FILENAME_LEN];
	snprintf(id, MAX_FILENAME_LEN, "%d", myid);

	// Initialize writing of logfile
	char logfname[MAX_FILENAME_LEN] = "\0";
	strAppend(logfname, dirname);
	strAppend(logfname, simname);
	strAppend(logfname, "logfile-node");
	strAppend(logfname, id);
	strAppend(logfname, ".txt");
	FILE *log;
	fileOpen(&log, logfname, "a");

	// Write log info
	char timeStampString[27] = "\0";
	getSystemTimeString(timeStampString);
	fprintf(log, "%s%s\n", timeStampString, str);
	fflush(log);
#if defined(__GNUC__)
	int fd = fileno(log);
	fsync(fd);
#elif defined(_MSC_VER)
	FlushFileBuffers(log);
#endif
	// Close log file
	fclose(log);
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
INLINE static int pickRndReaction(void) {
    probability	rate;
    int			i;
    probability	rnd, origRnd;
	static int	prevReact;

	// Choose reaction based on probability
	updateTree(prevReact);
	rate = state.reactProbTree[0];
    rnd = randomProb(2);

	rnd *= rate;
	origRnd = rnd;

	int curr = 1;
	probability curr_prob;
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
			printf("Fatal error: rate on node %d is %f. Exiting...\n\n", myid, rate);
			exit(EXIT_FAILURE);
		}
		else {
			return -1;
		}
	}

	if (DEBUGLEVEL >= 2) {
		if (i >= NO_OF_REACTIONS) {
			printf("pickRndReaction: ran out of reactions, i = %d no reactions = %d\n", i, NO_OF_REACTIONS);

			dumpReactProbTree();
			REACTION_PROBABILITY_TREE_INIT
				dumpReactProbTree();

			printf("\n original rnd no was %.40f\n", origRnd);
			exit(EXIT_FAILURE);
		}
	}
	
    return i;
}

/* Initialise the simulation
 */
static void initSysState(unsigned seed) {
	file_write_debug("Start initialization of the simulation");

	dsfmt_gv_init_gen_rand(seed);

	state.scaleFactor = 1.0;

    /* time, rate, particles, events, conversion, temperatures, etc. */
    state.time = 0.0;
	state.nextSynchTime = (ptime)(state.synchTime / 1000);
	state.events = 0;
	state.noMoreReactions = 0;
	state.noMoreReactionsLocal = 0;
	state.nextSynchEvents = (rcount)((float)state.synchEvents/(float)numprocs);
	state.localMonomerParticles = takeSome(GLOBAL_MONOMER_PARTICLES);
	state.conversion = 0;
	state.basetemp = BASETEMP;
	state.temp = STARTTEMP;
	state.deltatemp = state.temp - state.basetemp;

	/* moments of distribution */
	if (calcmoments) {
		state.momentDist[0] = 1;
		state.momentDist[1] = 1;
		state.momentDist[2] = 1;
	}
	
	/* free volume */
	if (calcfreevolume) {
		state.freeVolumeFraction = (VF0 + ALPHA_P * (state.temp - TG_P)) * state.conversion + (VF0 + ALPHA_M * (state.temp - TG_M)) * (1 - state.conversion);
	}

    /* MWD_INITS initialises those ms_cnts not set to 0. rely on malloc/memset for the rest */
    MWD_INITS

	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		state.ms_cnts[i] = takeSome(state.ms_cnts[i]);
	}
	
    /* initialise mwds */
    /* enter all mols into mwds */
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		for (int a = 0; a < MAX_ARMS; a++) {
	        state.mwds[i][a].maxEntries = START_MWD_SIZE;
	        state.mwds[i][a].mwd_tree = malloc((2 * START_MWD_SIZE - 1) * sizeof(pcount));
			memset((&state.mwds[i][a])->mwd_tree, 0, (2 * START_MWD_SIZE - 1) * sizeof(pcount));
		}
    }

	state.initialMonomerMolecules = monomerCount();
	state.currentMonomerMolecules = state.initialMonomerMolecules;
	state.volume = (state.localMonomerParticles / AVOGADRO) / MONOMERCONCENTRATION;

	REACTIONS_INIT

	int tmp[NO_OF_MOLSPECS] = ARMS_INIT;
	memcpy(state.arms, tmp, NO_OF_MOLSPECS * sizeof(int));
	memset(state.reactProbTree, 0, (2 * REACT_PROB_TREE_LEAVES - 1) * sizeof(probability));
	REACTION_PROBABILITY_TREE_INIT

	if (explicit) {
		state.expMols = malloc(sizeof(MoleculeList) * NO_OF_MOLSPECS);
		for (int i = 0; i < NO_OF_MOLSPECS; i++) {
			size_t initialSize = 10485760;
			if (i >= MAXSIMPLE) {
				state.expMols[i].mols = malloc(initialSize * sizeof(chainLen) * (size_t)state.arms[i]);
				state.expMols[i].maxMolecules = initialSize;
			}
		}
	}

	file_write_debug("Initialization of the simulation complete");
}

INLINE static int pickRndChainLen(pcount *mwd_tree, int size) {

	probability     prob, curr_prob;

	int             length;
	int             curr = 1;

	prob = randomProb(3);  
	prob *= mwd_tree[0];

	mwd_tree[0]--;
	while (curr < size) {
		curr = curr << 1;
		curr_prob = (probability)mwd_tree[curr - 1];
		if (prob >= curr_prob) {
			curr++;
			prob -= curr_prob;
		}
		if (DEBUGLEVEL >= 2) {
			if (mwd_tree[curr - 1] < 0) {
				return (-1);
			}
		}
		mwd_tree[curr - 1]--;
	}

	length = curr - size + 1;
	return (length);

}

/*
 *  Pick a random molecule which matches the spec and delete it from
 *  system.
 */
static int pickRndMolecule(const int spec_index, chainLen *lens, int *arms) {

	if (DEBUGLEVEL >= 2) {
		assert(state.ms_cnts[spec_index] > 0);
		if (state.ms_cnts[spec_index] <= 0) {
			printf("Trying to pick %s when there are none\n", name(spec_index));
			print_state();
			dumpReactProbTree();
			return -1;
		}
	}

    /* decrement global and spec local count */
	state.ms_cnts[spec_index]--;

    /* simple molecule, we're done */
    if (spec_index < MAXSIMPLE) {
        return 0;
    }
	else if (spec_index < MAXPOLY) {
		/* Poly species: pick one chain lengths according to MWD */
		*arms = 1;
		pcount *mwd_tree = state.mwds[spec_index][0].mwd_tree;
		int size = state.mwds[spec_index][0].maxEntries;
		lens[0] = pickRndChainLen(mwd_tree, size);

		if (calcmoments) {
			// Adjust moments of distribution for disappearance reactant
			state.momentDist[0]--;
			state.momentDist[1] -= lens[0];
			state.momentDist[2] -= lens[0] * lens[0];
		}

		if (DEBUGLEVEL >= 2) {
			if (lens[0] < 0) {
				printf("error: chain length < 0 for species %s(arm=0)\n", name(spec_index));
				print_state();
				dumpReactProbTree();
				abort();
			}
		}
	}
	else {
		/* Complex species: pick several chain lengths according to MWD */
		*arms = state.arms[spec_index];
		for (int a = 0; a < *arms; a++) {
			pcount *mwd_tree = state.mwds[spec_index][a].mwd_tree;
			int size = state.mwds[spec_index][a].maxEntries;
			lens[a] = pickRndChainLen(mwd_tree, size);

			if (DEBUGLEVEL >= 2) {
				if (lens[a] < 0) {
					printf("error: chain length < 0 for species %s(arm=%d)\n", name(spec_index), a);
					print_state();
					dumpReactProbTree();
					abort();
				}
			}
		}
		if (calcmoments) {
			// Adjust moments of distribution for disappearance reactant
			state.momentDist[0]--;
			chainLen comb_len = 0;
			for (int a = 0; a < *arms; a++) {
				comb_len += lens[a];
			}
			state.momentDist[1] -= comb_len;
			state.momentDist[2] -= comb_len * comb_len;
		}
	}

	return 0;
}

static int pickRndMolecule_order1(const int spec_index, chainLen *lens, int *arms) {

	if (DEBUGLEVEL >= 2) {
		if (state.ms_cnts[spec_index] <= 0) {
			printf("Trying to pick %s when there are none\n", name(spec_index));
			print_state();
			dumpReactProbTree();
			return -1;
		}
	}

    /* simple molecule, we're done */
    if (spec_index < MAXSIMPLE) {
    	state.ms_cnts[spec_index]--;
        return 0;
    }
	else {
		*arms = state.arms[spec_index];
		pcount lastIx = state.ms_cnts[spec_index]-1;
		pcount rndIx = (pcount)(randomProb(3) * (double)lastIx + 0.5); // Fast rounding by adding 0.5, then casting to int
		chainLen *pickedPos = &state.expMols[spec_index].mols[rndIx*(*arms)];
		size_t memSize = (*arms) * sizeof(chainLen);
		memcpy(lens, pickedPos, memSize);

		if (rndIx != lastIx) {
			memcpy(pickedPos, &state.expMols[spec_index].mols[lastIx*(*arms)], memSize);
		}

	    /* decrement global and spec local count */
	    state.ms_cnts[spec_index]--;

		if (calcmoments) {
			// Adjust moments of distribution for disappearance reactant
			state.momentDist[0]--;
			if (spec_index < MAXPOLY) {
				state.momentDist[1] -= (*lens);
				state.momentDist[2] -= (*lens) * (*lens);
			}
			else {
				chainLen comb_len = 0;
				for (int a = 0; a < (*arms); a++) {
					comb_len += lens[a];
				}
				state.momentDist[1] -= comb_len;
				state.momentDist[2] -= comb_len * comb_len;
			}
		}
		return 0;
	}
}

INLINE static const char *name(int index) {
    switch (index) 
        MOLECULENAMES
    
    return 0;
}


INLINE static const char *rname(int index) {
	switch (index) 
		REACTIONNAMES
	
	return 0;
}

/* Scales number of particles and (rate coefficients accordingly) by
 * factor.
 */
static void scaleSystem(double factor) {

	if (factor < 1.0) {
		printf("Scaling factors < 1 not yet dealt with\n");
		abort();
	}
	
	print_state();
	state.scaleFactor = factor;

	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		if (i < MAXSIMPLE) {
			state.ms_cnts[i] = (pcount)((double)state.ms_cnts[i] * factor);
			continue;
		}

		for (int a = 0; a < state.arms[i]; a++) {
			int j;
			for (j = 0; j < state.mwds[i][a].maxEntries - 1; j++) {
				state.mwds[i][a].mwd_tree[j] = (pcount)((double)state.mwds[i][a].mwd_tree[j] * factor);
			}

			pcount total = 0;
			for (j = state.mwds[i][a].maxEntries - 1; j < state.mwds[i][a].maxEntries * 2 - 1; j++) {
				state.mwds[i][a].mwd_tree[j] = (pcount)((double)state.mwds[i][a].mwd_tree[j] * factor);
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

	state.initialMonomerMolecules = (pcount)(state.initialMonomerMolecules * factor);
	state.currentMonomerMolecules = (pcount)(state.currentMonomerMolecules * factor);
	state.volume *= factor;

	printf("System scaled by factor of %f\n",factor);
	print_state();

}

/* Print reactants and results of a reaction
 *   (we could print the name of the reaction as well, would need to
 *   generate code for this -- see name (..))
 */
static void print_reaction(int reaction_index) {
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
static void react(void) {
	file_write_debug("  Start react() function; performing reactions until next sycnchronization point is reached");


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
		if (explicit) {
			prm1 = pickRndMolecule_order1(react1_ind, react1_lens, &react1_arms);
		}
		else {
			prm1 = pickRndMolecule(react1_ind, react1_lens, &react1_arms);
		}

		if (state.reactions[reactionIndex].arg_ms2 != NO_MOL) {
			react2_ind = state.reactions[reactionIndex].arg_ms2;
			if (explicit) {
				prm2 = pickRndMolecule_order1(react2_ind, react2_lens, &react2_arms);
			}
			else {
				prm2 = pickRndMolecule(react2_ind, react2_lens, &react2_arms);
			}

		}

		if (DEBUGLEVEL >= 2) {
			if (prm1 != 0 || prm2 != 0) {
				printf("Problem reaction %d\n", reactionIndex);
				abort();
			}
		}

		/*
		 * Given a reaction index, the lengths l1 and l2  of the reactants (0
		 * in case of Simple), determine type and length(s) of resulting
		 * molecule(s)
		 */
		switch (reactionIndex)
			DO_REACT_BODY

		if (DEBUGLEVEL >= 2) {
			assert(prod1_ind != -1);
		}

		state.ms_cnts[prod1_ind]++;
		if (prod1_ind >= MAXSIMPLE) {
			if (explicit) {
				adjustMolCnt_order1(prod1_ind, prod1_lens, prod1_arms);
			}
			else {
				adjustMolCnt(prod1_ind, prod1_lens, prod1_arms, 1);
			}

			if (DEBUGLEVEL >= 2) {
				assert(prod1_arms == state.arms[prod1_ind]);
			}
		}
		if (no_of_res == 2) {
			state.ms_cnts[prod2_ind]++;
			if (prod2_ind >= MAXSIMPLE) {
				if (explicit) {
					adjustMolCnt_order1(prod2_ind, prod2_lens, prod2_arms);
				}
				else {
					adjustMolCnt(prod2_ind, prod2_lens, prod2_arms, 1);
				}
			}
		}

		// Integer overflow protection (currently detects overflow on propagation, not combination)
		if ((prod1_lens[0] == CHAINLENLIMIT) ||
			(prod2_lens[0] == CHAINLENLIMIT)) {
			printf("\nError: a particle of species %s on node %d has reached its maximum chain length.\nWriting data and exiting...\n\n", name(prod1_ind), myid);
			printf("\nNode %d exited with errors\n", myid);
			exit(EXIT_FAILURE);
		}

		if (recalcconversion) {
			state.currentMonomerMolecules = monomerCount();
			state.conversion = conversion();
		}

		if (simulateheating) {
			state.deltatemp += (state.reactions[reactionIndex].energy / state.volume);
			state.temp = state.basetemp + state.deltatemp;
		}

		probability     rate;
		ptime			deltatime;

		rate = state.reactProbTree[0];

		deltatime = (-logRndArray()) / rate;
		state.time += deltatime;

		if (COOLINGRATE != 0) {
			state.deltatemp *= exp(-1 * coolingrate * deltatime);
			state.temp = state.basetemp + state.deltatemp;
		}

		if (calcfreevolume) {
			state.freeVolumeFraction = (VF0 + ALPHA_P * (state.temp - TG_P)) * state.conversion + (VF0 + ALPHA_M * (state.temp - TG_M)) * (1 - state.conversion);
		}

		if (DEBUGLEVEL >= 2) {
			assert(rate >= 0);
		}

		state.events++;
	} while (((state.time < state.nextSynchTime) || (state.synchTime == 0)) 
		    &&
		    ((state.events < state.nextSynchEvents) || (state.synchEvents == 0)));

    RANK printf("done!\n");

	file_write_debug("  React() has finished");
}

// Prints the mwds and molecule counts
static void print_state() {
	printf("\n\n");
    for (int i = 0; i < NO_OF_MOLSPECS; i++) {
        if (state.ms_cnts[i] > 0) {
            printf("Species %s (%llu): ", name(i), state.ms_cnts[i]);
            if (i < MAXSIMPLE) {
                printf("%llu\n", state.ms_cnts[i]);
            } else {
				pcount cnts[MAX_ARMS];
				int offset = state.mwds[i][0].maxEntries;
				for (int j = 0; j < offset; j++) {
					pcount total = 0;
					for (int a = 0; a<state.arms[i]; a++) {
						cnts[a] = getMolCnt(i, a, j);
						total += cnts[a];
					}
					if (total > 0) {
						printf("%d", j + 1);
						for (int a=0; a<state.arms[i]; a++) {
							printf("\t%llu", cnts[a]);
						}
						printf("\n");
					}
				}
            }
        }
    }
	printf("\n");
}

static void print_state_summary(const int m, const ptime *simtimes, const float *simconversions, const double *simtemps, const pcount *statecnts, const pcount *statemomentdists) {

	size_t nameLens[NO_OF_MOLSPECS];
	size_t maxNumLen = strlen("0.0000e+00");
	size_t nodeIDLen = strlen("Node XYZ ");
	size_t sum = 0;
	int nodesToPrint = 1;
	int rankLen;
	int distNo = 3;

	if (NODEVERBOSELEVEL >= 1) {
		nodesToPrint = numprocs;
	}

	sum += nodeIDLen;
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		nameLens[i] = max_value(strlen(name(i)), maxNumLen) + 1;
		sum += nameLens[i];
	}

	printf("\n");
	for (size_t i = 0; i < sum; i++) {
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
	printf("Wall time (s) = %f\n", readTimerSec(&state.wallTime));
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

	if (simulateheating) {
		// Temperature
		if (m == PRESTIRR) {
			for (int i = 0; i < nodesToPrint; i++) {
				printf("Temperature (K) on node %d = %.2f\n", i, (round(simtemps[i] * 100) / 100));
			}
		}
		else if (m == POSTSTIRR) {
			printf("Temperature (K) = %.2f\n", (round(simtemps[0] * 100) / 100));
		}
		printf("\n");
	}

	// Concentrations (and optionally moments of distribution, NACL, WACL, and PDI)
	for (int i = 0; i < nodeIDLen; i++)
		printf(" ");
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		printf("%s", name(i));
		size_t len_name = strlen(name(i));
		for (int j = 0; j < (nameLens[i] - len_name); j++) {
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
			size_t len_tmp = strlen(tmp);
			for (int k = 0; k < (nameLens[j] - len_tmp); k++)
				printf(" ");
		}
		printf("\n");
	}

	if (calcmoments) {
		printf("\n");
		for (size_t i = 0; i < nodeIDLen; i++)
			printf(" ");
		printf("0th moment 1st moment 2nd moment\n");
		for (size_t i = 0; i < nodeIDLen; i++)
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
	}

	for (size_t j = 0; j < sum; j++)
		printf("-");
	printf("\n");
	printf("\n");

}

static void state_summary(const int m) {

	int distNo = 3; // zeroth, first, and second moment of distribution

	// Collect simulation times
	ptime *workerStateTime = NULL;
	RANK workerStateTime = malloc(sizeof(ptime) * (size_t)numprocs);
	MPI_Gather(&state.time, 1, MPI_DOUBLE, workerStateTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Collect conversions
	float *workerStateConversion = NULL;
	RANK workerStateConversion = malloc(sizeof(float) * (size_t)numprocs);
	MPI_Gather(&state.conversion, 1, MPI_FLOAT, workerStateConversion, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

	// Collect simulation temperatures
	double *workerStateTemp = NULL;
	if (simulateheating) {
		RANK workerStateTemp = malloc(sizeof(double) * (size_t)numprocs);
		MPI_Gather(&state.temp, 1, MPI_DOUBLE, workerStateTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	// Collect species counts
	pcount *workerStateCnts = NULL;
	RANK workerStateCnts = malloc(sizeof(pcount) * (size_t)numprocs * NO_OF_MOLSPECS);
	MPI_Gather(&state.ms_cnts, NO_OF_MOLSPECS, MPI_UNSIGNED_LONG_LONG, workerStateCnts, NO_OF_MOLSPECS, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

	// Collect moments of distribution
	pcount *workerStateMomentDist = NULL;
	if (calcmoments) {
		RANK workerStateMomentDist = malloc(sizeof(pcount) * (size_t)numprocs * distNo);
		MPI_Gather(&state.momentDist, distNo, MPI_UNSIGNED_LONG_LONG, workerStateMomentDist, distNo, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	}

	// Print the state, then free up memory
	RANK print_state_summary(m, workerStateTime, workerStateConversion, workerStateTemp, workerStateCnts, workerStateMomentDist);

	free(workerStateTime);
	free(workerStateTemp);
	free(workerStateCnts);
	free(workerStateMomentDist);
}


/* Prints maxChainLen for all poly species */
static void printMaxChainLens(const int mode, const unsigned *workerMaxChainLens) {

	size_t maxNameLen = 10; // 9 characters and 1 space
	size_t nodeIDLen = 11; // length unsigned int (10 digits max) and 1 space
	size_t sum = 0;
	size_t nodesToPrint = 1;
	size_t rankLen = 1;
	size_t strLen = 1;

	if (NODEVERBOSELEVEL >= 1) {
		nodesToPrint = numprocs;
	}

	sum += maxNameLen;
	sum += nodeIDLen * nodesToPrint;

	printf("\n");
	for (size_t i = 0; i < sum; i++) {
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

	size_t pos = MAXSIMPLE;
	for (int i = MAXSIMPLE; i < MAXPOLY; i++) {
		//Print name + some spaces
		printf("%s", name(i));
		size_t max = (maxNameLen - strlen(name(i)));
		for (size_t j = 0; j < max; j++) {
			printf(" ");
		}

		//Print for each node
		for (size_t j = 0; j < nodesToPrint; j++) {

			// Print the value
			unsigned len = workerMaxChainLens[pos + j * MAXPOLY];
			printf("%u", len);

			// Print some more spaces
			if (len == 0) {
				strLen = 1;
			}
			else {
				strLen = (int)log10(len) + 1;
			}
			for (size_t k = 0; k < (nodeIDLen - strLen); k++) {
				printf(" ");
			}
		}
		printf("\n");
		pos++;
	}

	for (size_t i = 0; i < sum; i++) {
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
static size_t communicateMaxChainLens(const int mode, const unsigned *maxChainLens) {
	if (explicit) {  // Currently only compatible with explicit on

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

		if (VERBOSELEVEL >= 1) {
			// Collect maxChainLens from all nodes (poly species only), and print them
			unsigned *workerMaxChainLens = NULL;
			RANK workerMaxChainLens = malloc(sizeof(unsigned) * numprocs * MAXPOLY);
			MPI_Gather(maxChainLens, MAXPOLY, MPI_UNSIGNED, workerMaxChainLens, MAXPOLY, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

			RANK printMaxChainLens(mode, workerMaxChainLens);

			free(workerMaxChainLens);
		}

		if (VERBOSELEVEL >= 1) {
			RANK printf("communicateMaxChainLens: bytesNeeded = %zu\n", bytesNeeded);
		}
		return bytesNeeded;
	}
	else {
		return 0;
	}
}

static void print_kinetic_model(void) {
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
	for (int i = 0; i < NO_OF_REACTIONS; i++) {
		reaction *p = state.reactions + i;
		char *p1 = (char*)((p->res_ms1!=NO_MOL)?name(p->res_ms1):" ");
		char *p2 = (char*)((p->res_ms2!=NO_MOL)?name(p->res_ms2):" ");
		char *r1 = (char*)((p->arg_ms1!=NO_MOL)?name(p->arg_ms1):" ");
		char *r2 = (char*)((p->arg_ms2!=NO_MOL)?name(p->arg_ms2):" ");
		double k;
		if (p->arg_ms2 == NO_MOL) { // Unimolecular reaction
			k = state.reactions[i].rc;
		}
		else if (p->arg_ms1 == p->arg_ms2) { // Bimolecular reaction with identical reactants
			k = state.reactions[i].rc * state.volume * AVOGADRO / SZYMANSKI;
		}
		else { // Bimolecular reaction with nonidentical reactants
			k = state.reactions[i].rc * state.volume * AVOGADRO;
		}
		printf("%d:\t%s\t+\t%s\t-->\t%s\t+\t%s\t k = %.5e\n", i, r1, r2, p1, p2, k);
	}
}

static pcount getConvertedMonomer(void) {
    int i, j;
	pcount convertedMonomer = 0;

	if (explicit) {
		for (i = 0; i < NO_OF_MOLSPECS; i++) {
			if (i >= MAXSIMPLE) {
				for (j = 0; j < state.arms[i] * state.ms_cnts[i]; j++) {
					convertedMonomer += state.expMols[i].mols[j];
				}
			}
		}
	}
	else {
		int length, offset;
		for (i = 0; i < NO_OF_MOLSPECS; i++) {
			if (i >= MAXSIMPLE) {
				for (int a = 0; a < state.arms[i]; a++) {
					length = 0;
					offset = state.mwds[i][a].maxEntries - 1;
					for (j = offset; j < 2 * offset + 1; j++) {
						length++;
						convertedMonomer += length*state.mwds[i][a].mwd_tree[j];
					}
				}
			}
		}
	}

	return convertedMonomer;
}

static void monomerAuditLocal(const char *str) {
	if (!monomeraudit) {
		return; // Do not perform
	}

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
				printf("\nLocal monomer audit (%s) failed on node %d! Discrepancy = %llu\n\n", str, i, workerDiscrepancies[i]);
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

/* Performs a+b=c for an array
 */
INLINE static void dvecAdd(pcount *c, pcount *a, pcount *b, unsigned n) {
	for (size_t i = 0; i < n; i++) {
			c[i] = a[i] + b[i];
	}
}

INLINE static int sizeExceeded(StatePacket *sp) {
    return(sp->stateTooBig);
}


static void buildTreeFromLeaves(pcount *tree, pcount *leaves, int maxLeaves) {
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
static size_t getMaxChainLens(unsigned *max_len_arr_ptr) {
	size_t bytesNeeded = 0, pos = 0;
	
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		if (i < MAXSIMPLE) {
			max_len_arr_ptr[pos] = 0;
			pos++;
		} else {
			if (explicit) {
				if (state.arms[i] == 1) { // For single arm species, store max chain length as usual
					unsigned maxLen = 0;
					for (pcount x = 0; x < state.ms_cnts[i]; x++) {
						unsigned tmp = state.expMols[i].mols[x];
						if (tmp > maxLen)
							maxLen = tmp;
					}
					max_len_arr_ptr[pos] = maxLen;
					if (VERBOSELEVEL >= 3) {
						RANK printf("max_len_arr_ptr[%s] on node %d = %u\n", name(i), myid, maxLen);
					}
					bytesNeeded += max_len_arr_ptr[pos];
					pos++;
				}
				else { // For complex species, 1st entry stores number of molecules, rest are zero
					max_len_arr_ptr[pos] = (unsigned)state.ms_cnts[i];
					if (VERBOSELEVEL >= 3) {
						RANK printf("max_len_arr_ptr[%s] on node %d = %u\n", name(i), myid, max_len_arr_ptr[pos]);
					}
					pos++;
					for (int a = 1; a < state.arms[i]; a++) {
						max_len_arr_ptr[pos++] = 0;
					}
				}
			}
			else {
				int currMax;
				for (int a = 0; a < state.arms[i]; a++) {
					int offset = state.mwds[i][a].maxEntries - 1;
					currMax = offset - 1;
					for (int j = offset; j < 2 * offset + 1; j++) {
						if (state.mwds[i][a].mwd_tree[j] > 0) {
							currMax = j;
						}
					}
					max_len_arr_ptr[pos] = currMax - offset + 1;
					bytesNeeded += max_len_arr_ptr[pos];

					pos++;
				}
			}
		}
    }
	bytesNeeded *= sizeof(pcount); // Calculate space needed to store MWDs for all species
	if (VERBOSELEVEL >= 1) {
		RANK printf("getMaxChainLens on node %d: bytesNeeded = %zu\n", myid, bytesNeeded);
	}
	return bytesNeeded;
}

/*
 * Calculate number of bytes required for the header of the state packet
 */
static size_t stateCommHeaderBytes(void) {

	return (sizeof(StatePacket)
		  + sizeof(pcount)*NO_OF_MOLSPECS // Mol counts
		  + sizeof(rcount)*NO_OF_REACTIONS // Reaction event counts
          + sizeof(unsigned)*TOTAL_ARMS);  // Maximum chain lengths
}

static size_t requiredStateCommSize(void) {
	// Allocate memory, get maxChainLens
	unsigned *maxChainLens = malloc(sizeof(unsigned) * TOTAL_ARMS);
	getMaxChainLens(maxChainLens);
	
	// communicate MaxChainLens between workers, finds largest numbers, calculate space needed for storing MWDs
	size_t bytesForMWDs = communicateMaxChainLens(PRESTIRR, maxChainLens);

	free(maxChainLens);
	return stateCommHeaderBytes() + bytesForMWDs;
}

static void stateToComm(StatePacket **outStatePacket, StatePacket **inStatePacket) {
	file_write_debug("    Starting stateToComm() function");

	// Packet format
	// [ Header                                                  | Body                                 ]
	// [ StatePacket    | pcount[]  | rcount[]    | int[]        | pcount[] | pcount[] | ... | pcount[] ]
	// [ outStatePacket | molCounts | reactCounts | maxChainLens | packetMWDs ...                    E  ]
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

	// Define start of array containing total number of reaction events, then copy event counts
	pcount *reactCounts = (rcount*)(molCounts + NO_OF_MOLSPECS);
	memcpy(reactCounts, state.react_cnts, sizeof(rcount) * NO_OF_REACTIONS);

	// Define start of array containing maximum chain lengths for each species, then copy maximum chain lengths
	unsigned *maxChainLens = (unsigned*)(reactCounts + NO_OF_REACTIONS);
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
			if (explicit) {
				// First set count for all chainlengths to 0, then populate the MWD packet for this species
				if (state.arms[i] == 1) {
					memset(packetMwds, 0, maxChainLens[pos] * sizeof(pcount));
					for (int m = 0; m < state.ms_cnts[i]; m++) {
						chainLen len = state.expMols[i].mols[m];
						packetMwds[len - 1]++;
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
					RANK if (VERBOSELEVEL >= 3) {
						printf("Before stirr: %llu molecules of %s in packet on node 0", state.ms_cnts[i], name(i));
						if (state.ms_cnts[i] > 0) {
							printf("; arm lengths:");
							for (int t = 0; t < state.ms_cnts[i] * state.arms[i]; t++) {
								printf(" %u", *((chainLen*)packetMwds + t));
							}
						}
						printf("\n");
					}
					packetMwds = (pcount*)((char*)packetMwds + bytes);
					pos++;

					for (int a = 1; a < state.arms[i]; a++) {
						maxChainLens[pos++] = 0;
					}
				}
			}
			else {
				for (int a = 0; a < state.arms[i]; a++) {
					pcount *leaves = state.mwds[i][a].mwd_tree + state.mwds[i][a].maxEntries - 1;
					memcpy(packetMwds, leaves, maxChainLens[pos] * sizeof(pcount));
					packetMwds += maxChainLens[pos];
					pos++;
				}
			}
		}
	}

	file_write_debug("    stateToComm() has finished");
}

static int comparitorComplex(const void *m1_, const void *m2_) {
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
static void commToState(StatePacket **inStatePacket) {
	file_write_debug("    Starting commToState() function");

	RANK printf("commToState running...\n");
	
	state.noMoreReactions = (*inStatePacket)->noMoreReactions;
    state.time = (*inStatePacket)->time / numprocs;
	state.deltatemp = (*inStatePacket)->deltatemp / numprocs;

	pcount *speciesCounts = (pcount*)((*inStatePacket) + 1);
	pcount *reactCounts = (rcount*)(speciesCounts + NO_OF_MOLSPECS);
	unsigned *maxChainLens = (unsigned*)(reactCounts + NO_OF_REACTIONS);
	pcount *mwd = (pcount*)(maxChainLens + TOTAL_ARMS);

	// Send count of each reaction event on node 0, set counts on other nodes to 0
	for (int i = 0; i < NO_OF_REACTIONS; i++) {
		RANK {
			state.react_cnts[i] = reactCounts[i];
		}
		else {
			state.react_cnts[i] = 0;
		}
	}

	// Copy mwd trees into state and take some particles
	int pos = 0;
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {
		if (i < MAXSIMPLE) {
			// divide particles up amongst processes
			state.ms_cnts[i] = takeSome(speciesCounts[i]);
			pos++;
			continue;
		}

		state.ms_cnts[i] = 0;
		
		if (explicit) {
			if (state.arms[i] == 1) { // poly species

				unsigned maxChainLen = maxChainLens[pos];
				if (VERBOSELEVEL >= 2) {
					RANK printf("maxChainLen[%s] = %u\n", name(i), maxChainLen);
				}
				int takeInterval = numprocs;
				int takeCounter = myid; // myid functions as offset
				pcount molsAddedSoFar = 0;

				for (unsigned j = 0; j < maxChainLen; j++) {
					pcount mols = *(mwd++);
					chainLen length = j + 1;
					while (mols > 0) {
						if (molsAddedSoFar >= state.expMols[i].maxMolecules) {
							RANK if (VERBOSELEVEL >= 3) {
								printf("masf = %llu mm = %llu\n", molsAddedSoFar, state.expMols[i].maxMolecules);
							}
							printf("system state expanded for %s in commToState\n", name(i));
							state.expMols[i].maxMolecules = (pcount)((double)state.expMols[i].maxMolecules * 1.5);
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
				if (VERBOSELEVEL >= 3) {
					RANK printf("After stirr: %llu molecules of %s on node 0\n", molecules, name(i));
				}

				pcount end = start + molecules;
				chainLen *list = (chainLen*)mwd;

				// The complex species must be sorted! This is because different processes
				// will merge them into the packet in different orders resulting in problems
				// when the different processes then take molecules to resume simulation with.
				currentComparisonComplexity = state.arms[i];
				qsort(list, maxChainLens[pos], sizeof(chainLen)*(size_t)state.arms[i], comparitorComplex);

				state.ms_cnts[i] = 0;
				for (pcount c = start; c < end; c++) {
					chainLen *mol = list + state.arms[i] * c;
					state.ms_cnts[i]++;
					adjustMolCnt_order1(i, mol, state.arms[i]);
				}

				mwd = (pcount*)((char*)mwd + sizeof(chainLen) * maxChainLens[pos] * state.arms[i]);

				pos += state.arms[i];
				RANK if (VERBOSELEVEL >= 3) {
					printf("After stirr and sort: %llu molecules of %s on node 0", state.ms_cnts[i], name(i));
					if (state.ms_cnts[i] > 0) {
						printf("; arm lengths:");
						for (int t = 0; t < state.ms_cnts[i] * state.arms[i]; t++) {
							printf(" %u", (unsigned)state.expMols[i].mols[t]);
						}
					}
					printf("\n");
				}

			}
		}
		else {
			pcount ms_cnts_tmp_min = state.localMonomerParticles; // for some reason LONG_LONG_MAX is not defined
			pcount ms_cnts_tmp_max = 0;

			for (int a = 0; a < state.arms[i]; a++) {

				while (maxChainLens[pos] > state.mwds[i][a].maxEntries) {
					fprintf(stderr, "commToState: expanding tree memory for species %s(arm=%d)\n", name(i), a);
					free(state.mwds[i][a].mwd_tree);
					state.mwds[i][a].maxEntries *= 2;
					state.mwds[i][a].mwd_tree = malloc(sizeof(pcount) * (2 * state.mwds[i][a].maxEntries - 1));
				}
				memset(state.mwds[i][a].mwd_tree, '\0', (2 * state.mwds[i][a].maxEntries - 1) * sizeof(pcount));

				// new MWD becomes leaves in probability tree
				unsigned entries = state.mwds[i][a].maxEntries;
				pcount *leaves = state.mwds[i][a].mwd_tree + entries - 1;

				unsigned maxChainLen = maxChainLens[pos];
				assert(maxChainLen <= entries);

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

				buildTreeFromLeaves(state.mwds[i][a].mwd_tree, leaves, entries);
				pos++;
			}
			state.ms_cnts[i] = ms_cnts_tmp_min;
			if (ms_cnts_tmp_max - ms_cnts_tmp_min > 0)
				printf("Warning: %lld molecules of species %s temporarily lost upon communication\n", ms_cnts_tmp_max - ms_cnts_tmp_min, name(i));
		}
	}

	if (VERBOSELEVEL >= 1) {
		communicateMaxChainLens(POSTSTIRR, maxChainLens);
	}

	REACTION_PROBABILITY_TREE_INIT
	
	// Global monomer audit
	if (monomeraudit) {
		RANK { 
			if ((*inStatePacket)->globalAllMonomer != GLOBAL_MONOMER_PARTICLES) {
				printf("\nGlobal monomer audit has failed! Discrepancy = %llu\n\n", (*inStatePacket)->globalAllMonomer - GLOBAL_MONOMER_PARTICLES);
			}
			else {
				printf("Global monomer audit passed.\n");
			}
		}
		// setup data for local monomer audits during the next simulation phase
		state.currentMonomerMolecules = monomerCount();
		state.initialMonomerMolecules = state.currentMonomerMolecules + getConvertedMonomer();
		if (VERBOSELEVEL >= 3) {
			RANK printf("initialMonomerMolecules = %lld\n", state.initialMonomerMolecules);
		}
	}

	file_write_debug("    commToState() has finished");
}

static void checkMWDs(pcount *mwd, unsigned totalLength, char *str) {
	for (unsigned i = 0; i < totalLength; i++) {
		if (mwd[i] > state.localMonomerParticles)
			printf("******************* %s: Bad data after stirring totalLen = %u\n", str, totalLength);
	}
}

/* Combines the state information which has been sent. Cyber 
 * equivalent of stirring.
 */
static void stirr(void *in_, void *inout_, int *len, MPI_Datatype *datatype) {
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
	pcount *reactCounts_in = (rcount*)(specCounts_in + NO_OF_MOLSPECS);
	unsigned *maxChainLens_in = (unsigned*)(reactCounts_in + NO_OF_REACTIONS);
	pcount *mwd_in = (pcount*)(maxChainLens_in + TOTAL_ARMS);

	// Define layout of packet 'inout'
	pcount *specCounts_inout = (pcount*)(inout+1);
	pcount *reactCounts_inout = (rcount*)(specCounts_inout + NO_OF_MOLSPECS);
	unsigned *maxChainLens_inout = (unsigned*)(reactCounts_inout + NO_OF_REACTIONS);
	pcount *mwd_inout = (pcount*)(maxChainLens_inout + TOTAL_ARMS);

	pcount *mwd_inout_start = mwd_inout;
	pcount *mwd_tmp = malloc(stateCommSize);
	pcount *mwd_pos = mwd_tmp;
	size_t totalLength = 0, miscBytes = 0;

	// Merge reaction event counts
	for (int i = 0; i < NO_OF_REACTIONS; i++) {
		reactCounts_inout[i] += reactCounts_in[i];
	}

	// Merge species counts and MWDs
	int pos = 0;
	for (int i = 0; i < NO_OF_MOLSPECS; i++) {

		specCounts_inout[i] += specCounts_in[i];
		if (i < MAXSIMPLE) {
			maxChainLens_inout[i] = 0;
			pos++;
			continue;
		}

		if (explicit) {
			if (i >= MAXPOLY) { // Complex species
				chainLen *mwd_pos_old = (chainLen*)mwd_pos;

				size_t inoutBytes = maxChainLens_inout[pos] * sizeof(chainLen) * state.arms[i];
				memcpy(mwd_pos, mwd_inout, inoutBytes);
				mwd_pos = (pcount*)((char*)mwd_pos + inoutBytes);

				size_t inBytes = maxChainLens_in[pos] * sizeof(chainLen) * state.arms[i];
				memcpy(mwd_pos, mwd_in, inBytes);
				mwd_pos = (pcount*)((char*)mwd_pos + inBytes);

				maxChainLens_inout[pos] += maxChainLens_in[pos];

				if (VERBOSELEVEL >= 1) {
					RANK{
						printf("stirr: %u molecules of %s in stirred packet", maxChainLens_inout[pos], name(i));
						if (maxChainLens_inout[pos] > 0) {
							printf("; arm lengths:");
							for (unsigned t = 0; t < maxChainLens_inout[pos] * state.arms[i]; t++) {
								printf(" %u", (unsigned)mwd_pos_old[t]);
							}
						}
						printf("\n");
					}
				}

				miscBytes += inoutBytes + inBytes;
				pos += state.arms[i];
				mwd_inout = (pcount*)((char*)mwd_inout + inoutBytes);
				mwd_in = (pcount*)((char*)mwd_in + inBytes);
			}
			else { // Poly species
				/* Allow for the situation which will arise in which a mwd tree in some node
				   has undergone a size increase before the corresponding mwd in another node.

				   Note: this may no longer be possible now that nodes first agree on an adequate MWD tree size
				   Redundant code?
				 */
				unsigned common = (unsigned)min_value(maxChainLens_in[pos], maxChainLens_inout[pos]);
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
				checkMWDs(mwd_pos, maxChainLens_inout[pos], "mwd_pos");
				mwd_pos += maxChainLens_inout[pos];
				totalLength += maxChainLens_inout[pos];
				if (VERBOSELEVEL >= 2) {
					RANK printf("stirr: maxChainLens_inout[pos] for species %d is %u\n", i, maxChainLens_inout[pos]);
				}
				pos++;
			}
		}
		else {
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
		}
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

static int MPI_Allreduce_wrapper(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm ) {
	file_write_debug("    Start MPI reduction (stirring)");
	reduces++;
	Timer t;
    startTimer(&t);
    int r = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);

	unsigned long long rtime = readTimer(&t);
	total_rtime += rtime;
	
    RANK printf("Reduce time (us) = %llu for %zu bytes\n", rtime, stateCommSize);
	lastReduceTime = rtime;
	file_write_debug("    MPI reduction (stirring) complete");
    return r;
}

/* Dump state content -- currently not functional yet (data likely needs to be serialized) */
static void dumpStatePacket(StatePacket **inStatePacket, int stateCommSize) {

	// Open the file
	char fname[MAX_FILENAME_LEN] = "\0";
	strAppend(fname, dirname);
	strAppend(fname, simname);
	strAppend(fname, "stateDump");
	FILE *dump;
	fileOpen(&dump, fname, "w");

	// Check if file could be opened
	if (dump == NULL) {
		printf("Could not write to file %s\n", fname);
		exit(EXIT_FAILURE);
	}

	// Write to the file
	size_t r1 = fwrite(*inStatePacket, stateCommSize, 1, dump);
	printf("wrote %zu elements out of %d requested\n", r1, stateCommSize);

	//Close the file
	fclose(dump);
}

// ************* end communication code **********


// ************* start testing code *********

/* We're going to test a few things. */  

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
	char errormsg[100];
	snprintf(errormsg, 100, "PRNG_ARRAY_SIZE may not be smaller than %d", DSFMT_N64);
	mu_assert(PRNG_ARRAY_SIZE >= DSFMT_N64, errormsg);
	mu_assert((PRNG_ARRAY_SIZE % 2) == 0, "PRNG_ARRAY_SIZE must be an even number");

	// 1) Seed the PRNG, 2) test five consecutive values, 3) empty the PRNG arrays
	// Test path 0
	dsfmt_gv_init_gen_rand(1);
	mu_assert_double_eq(0.1193544251137, randomProb(0));
	mu_assert_double_eq(0.9124176151803, randomProb(0));
	mu_assert_double_eq(0.5031786702429, randomProb(0));
	mu_assert_double_eq(0.8712546575055, randomProb(0));
	mu_assert_double_eq(0.5324328025691, randomProb(0));
	randomProb(4);

	// Test path 2
	dsfmt_gv_init_gen_rand(1);
	mu_assert_double_eq(0.8806455748863, randomProb(2));
	mu_assert_double_eq(0.0875823848197, randomProb(2));
	mu_assert_double_eq(0.4968213297571, randomProb(2));
	mu_assert_double_eq(0.1287453424945, randomProb(2));
	mu_assert_double_eq(0.4675671974309, randomProb(2));
	randomProb(4);

	// Test path 3 -- currently same as path 0
	/*dsfmt_gv_init_gen_rand(1);
	mu_assert_double_eq(0.1193544251137, randomProb(0));
	mu_assert_double_eq(0.9124176151803, randomProb(0));
	mu_assert_double_eq(0.5031786702429, randomProb(0));
	mu_assert_double_eq(0.8712546575055, randomProb(0));
	mu_assert_double_eq(0.5324328025691, randomProb(0));
	randomProb(4);*/
}

/*  Test whether START_MWD_SIZE is a power of two */
MU_TEST(START_MWD_SIZE_test) {
	int value = START_MWD_SIZE;
	mu_assert(value > 0, "START_MWD_SIZE is not a power of two");
	int value_test = value & (~value + 1);
	mu_assert(value == value_test, "START_MWD_SIZE is not a power of two");
}

/* Tests to run */
MU_TEST_SUITE(test_suite) {

	MU_SUITE_CONFIGURE(&test_setup, &test_teardown);

	MU_RUN_TEST(myid_test);
	MU_RUN_TEST(dSFMT_test);
	MU_RUN_TEST(START_MWD_SIZE_test);
}

// ************* end testing code *********


/*  Parse and apply command line options */
static void argumentParsing(int argc, const char *argv[], unsigned *seed) {
	file_write_debug("Start parsing of command line arguments");

	unsigned arg_seed = UINT_MAX;
	unsigned arg_synchevents = UINT_MAX;
	unsigned arg_synchtime = UINT_MAX;
	char *arg_dirname = NULL;
	char *arg_simname = NULL;
	struct argparse_option options[] = {
		OPT_HELP(),
		OPT_GROUP("Basic options"),
		OPT_UINT('s', "seed", &arg_seed, "seed for PRNG (allowed range: 0 - 4294967294))"),
		OPT_UINT('e', "events", &arg_synchevents, "number of events between each two synchronizations (allowed range: 0 - 4294967294)"),
		OPT_UINT('t', "simtime", &arg_synchtime, "simulation time (ms) between each two synchronizations (allowed range: 0 - 4294967294)"),
		OPT_STRING('o', "output", &arg_dirname, "output directory"),
		OPT_STRING('n', "simname", &arg_simname, "name of the case (used as a prefix for filenames of output files"),
		OPT_END(),
	};
	struct argparse argparse;
	argparse_init(&argparse, options, usages, 0);
	argparse_describe(&argparse, "\nSimply - kinetic Monte Carlo simulator for polymerizations.", "\n");
	argc = argparse_parse(&argparse, argc, argv);
	if (arg_seed != UINT_MAX) {
		*seed = arg_seed + myid;
		RANK printf("Overriding compiled PRNG seed with user specified seed: %u\n", *seed);
	}
	else {
		RANK printf("PRNG seed: %u\n", *seed);
	}
	if (arg_synchevents != UINT_MAX) {
		state.synchEvents = (rcount)arg_synchevents;
		RANK printf("Overriding compiled number of events synchronization interval with user specified number: %llu\n", state.synchEvents);
	}
	else {
		RANK printf("Number of events between each two synchronizations: %llu\n", state.synchEvents);
	}
	if (arg_synchtime != UINT_MAX) {
		state.synchTime = (unsigned long long)arg_synchtime;
		RANK printf("Overriding compiled simulation time synchronization interval with user specified time (ms): %llu\n", state.synchTime);
	}
	else {
		RANK printf("Simulation time between each two synchronizations: %llu\n", state.synchTime);
	}
	if (arg_dirname != NULL) {
		parseDirname(arg_dirname);
		strAppend(dirname, arg_dirname);
		printf("User defined output path: %s\n", dirname);
	}
	if (arg_simname != NULL) {
		strAppend(simname, arg_simname);
		printf("User defined filename prefix: %s\n", simname);
		strAppend(simname, "-");
	}

	// Check if at least one sync interval is nonzero -- to be used in the future
	if ((state.synchTime == 0) && (state.synchEvents == 0)) {
		RANK printf("\nNo valid synchronization interval has been provided. Exiting...\n");
		exit(EXIT_FAILURE);
	}

	file_write_debug("Parsing of command line arguments complete");
}

static void compute(void) {
	file_write_debug("Starting simulation");

	RANK file_write_state(START);
	RANK file_write_state(PROFILES);

	StatePacket *outStatePacket = NULL;
	StatePacket *inStatePacket = NULL;

/* Stop if all of the following criteria are met:
    - max simtime is exceeeded
    - max number of events is exceeeded
    - max conversion number is exceeeded
   In addition, always stop when the maximum
   walltime has been exceeeded or no more reactions
   are possible on one or more nodes */
	while (((state.time       < MAX_SIM_TIME  )
		||  (state.events     < MAX_EVENTS    ) 
		||  (state.conversion < MAX_CONVERSION))
		&& (readTimerSec(&state.wallTime) < MAX_WALL_TIME * 60)
		&& (state.noMoreReactions == 0))
	{
		Timer work;
		startTimer(&work);
		rcount prev_events = state.events;

		react();

		unsigned long long wtime = readTimer(&work);
		total_wtime += wtime;
		RANK printf("Total calculation walltime (us) = %llu\n", total_wtime);
		state.currentMonomerMolecules = monomerCount();
		state.conversion = conversion();

		monomerAuditLocal("post reactions");

		RANK printf("Work time (us) on node %d = %llu\n", myid, wtime);
		RANK printf("Total number of events on node %d = %llu\n", myid, state.events);
		RANK printf("Computational speed on node %d = %.4f events/us\n", myid, (float)(state.events-prev_events)/(float)wtime);
		
		state_summary(PRESTIRR);

		RANK printf("About to synch, time = %f\n", state.time);

		file_write_debug("  Starting synchronization");

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
		free(outStatePacket);
		free(inStatePacket);

        // Create packet to send, and allocate memory for received packet
		outStatePacket = (StatePacket*)malloc(stateCommSize);
		inStatePacket = (StatePacket*)malloc(stateCommSize);

		stateToComm(&outStatePacket,&inStatePacket);

		/* Explain to MPI how type is defined.
		   This must be redone each time since each packet is a single
		   value whose size changes with each send. */
		MPI_Type_contiguous(stateCommSize, MPI_CHAR, &state_t); 
		MPI_Type_commit(&state_t); 
		MPI_Op_create(stirr, True, &myOp);
			
		// Perform reduction
		reduceRes = MPI_Allreduce_wrapper(outStatePacket, inStatePacket, 1, state_t, myOp, MPI_COMM_WORLD);
		if (reduceRes != MPI_SUCCESS) {
			RANK printf("\nError during synchronization! PMPI_Allreduce error code is %d\nExiting...\n\n", reduceRes);
			exit(EXIT_FAILURE);
		}

		// Free operation
		MPI_Op_free(&myOp);

		// Update the state after stirring
		commToState(&inStatePacket);
		file_write_debug("  Synchronization complete");
		monomerAuditLocal("post communicate");
	
		// Generate new seed
		if (changeseed) {
			unsigned nextSeed = (unsigned)(lastReduceTime * randomProb(0) + getpid() * randomProb(0));
			dsfmt_gv_init_gen_rand(nextSeed);
			RANK printf("New random seed on node 0 is %u\n", nextSeed);
		}

		// Recalculate what's needed
		state.temp = state.basetemp + state.deltatemp;
		state.conversion = conversion();
		if (calcfreevolume) {
			state.freeVolumeFraction = (VF0 + ALPHA_P * (state.temp - TG_P)) * state.conversion + (VF0 + ALPHA_M * (state.temp - TG_M)) * (1 - state.conversion);
		}
		if (calcmoments) {
			// CURRENTLY ONLY FUNCTIONAL WITH EXPLICIT STATE!
			state.momentDist[0] = 1;
			state.momentDist[1] = 1;
			state.momentDist[2] = 1;

			for (int i = 0; i < NO_OF_MOLSPECS; i++) {

				if (i >= MAXSIMPLE) {
					state.momentDist[0] += state.ms_cnts[i];
					size_t arms = state.arms[i];

					for (int j = 0; j < state.ms_cnts[i]; j++) {
						chainLen *lengths = &state.expMols[i].mols[arms*j];
						chainLen totalLen = 0;

						for (size_t a = 0; a < arms; a++) {
							totalLen += lengths[a];
						}
						state.momentDist[1] += totalLen;
						state.momentDist[2] += totalLen * totalLen;
					}
				}
			}
		}

		// Print synced results to screen and file
		state_summary(POSTSTIRR);
		RANK file_write_state(PROFILES);

		state.nextSynchTime += (ptime)(state.synchTime / 1000);
		state.nextSynchEvents += (rcount)((float)state.synchEvents/(float)numprocs);
    }
	file_write_debug("Simulation complete");

    return;
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#if defined(_MSC_VER)
	QueryPerformanceFrequency(&frequency);
#endif
	startTimer(&state.wallTime);

	retreiveHostname();
	RANK{
		printf("Starting up on host \"%s\"\n", hostname);
		printf("\n");
		printf("Simply version 0.99 beta prerelease\n");
		printf("Program compiled at %s on %s\n",__TIME__,__DATE__);
		printf("\n");
		printf("Number of nodes = %d\n", numprocs);
		printf("\n");
	}
	file_write_debug("");
	file_write_debug("Simply started");

	// Define arrays for PRNG
#if defined(__GNUC__)
	w128_t *dummy0 = memalign(16, (PRNG_ARRAY_SIZE / 2 + 1)*sizeof(w128_t));
	w128_t *dummy2 = memalign(16, (PRNG_ARRAY_SIZE / 2 + 1) * sizeof(w128_t));
	w128_t *logDummy2 = memalign(16, (PRNG_ARRAY_SIZE / 2 + 1) * sizeof(w128_t));
#elif defined (_MSC_VER) // MSVC automatically aligns _m128* types on 16-byte boundaries
	w128_t *dummy0 = malloc((PRNG_ARRAY_SIZE / 2 + 1) * sizeof(w128_t));
	w128_t *dummy2 = malloc((PRNG_ARRAY_SIZE / 2 + 1) * sizeof(w128_t));
	w128_t *logDummy2 = malloc((PRNG_ARRAY_SIZE / 2 + 1) * sizeof(w128_t));
#endif
	rndArray0 = (double *)dummy0;
	rndArray2 = (double *)dummy2;
	logRndArray2 = (double *)logDummy2;

	// Unit testing
	RANK {
		file_write_debug("Starting unit testing");
		MU_RUN_SUITE(test_suite);
		MU_REPORT();
		file_write_debug("Unit testing complete");
	}
	
	state.synchTime = SYNCH_TIME_INTERVAL;
	state.synchEvents = SYNCH_EVENTS_INTERVAL;
#ifndef SEED
    unsigned seed = getpid();
#else
	unsigned seed = SEED + myid;
#endif

	argumentParsing(argc, argv, &seed);

	// Initialize the calculation
	RANK printf("\n");
    initSysState(seed);
	RANK print_kinetic_model();
    RANK print_state();

	// Run the simulation
    compute();

	// Simulation is finished
	RANK {
		if (state.noMoreReactions != 0) {
			printf("Exiting prematurely as %d node(s) ran out of reactions to perform\n\n", state.noMoreReactions);
		}

		printf("Wall time (s) = %.2f\n", readTimerSec(&state.wallTime));
		printf("Computational time (us): chatting = %llu (avg = %llu), working = %llu\n", total_rtime, (reduces > 0 ? (total_rtime / reduces) : 0), total_wtime);
		printf("Averaged computational speed on node %d = %.4f events/us\n", myid, (float)(state.events) / (float)total_wtime);
		printf("Parallel efficiency = %.1f\n", (float)total_wtime / (float)(total_wtime + total_rtime) * 100);
		printf("\n");

		rcount totalEvents = 0;
		for (int i = 0; i < NO_OF_REACTIONS; i++) {
			totalEvents += state.react_cnts[i];
		}
		printf("Total number of events (on all nodes) = %llu\n", totalEvents);
		printf("Final conversion = %f\n", state.conversion);
		printf("Final simulation time = %f\n", state.time);
		printf("\n");
	}

	if (explicit) {
		compressState();
	}
	mergeMWDs();
	RANK file_write_MWDs();

	MPI_Finalize();

	exit(EXIT_SUCCESS);
}


