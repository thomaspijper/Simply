/* 
 *  Copyright (c) 2005 Gabriele Keller
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
 * 
 */

#include "genpolymer.h"
#include <inttypes.h>
#if defined(__GNUC__)
  #include <sys/time.h>
#endif

#define MAX_DATA_FILE 1000
#define RATES_VEC_SIZE 10
typedef double probability;
typedef double react_prob;
typedef double ptime;
typedef double rateCoefficient;
#ifdef LONGCHAINSUPPORT
typedef unsigned chainLen;
#define CHAINLENLIMIT UINT_MAX
#else
typedef short chainLen;
#define CHAINLENLIMIT SHRT_MAX
#endif
typedef unsigned long long pcount; // stores a number of particles

typedef struct {
	rateCoefficient rc;
	int arg_ms1;
	int arg_ms2;
	int res_ms1;
	int res_ms2;    
	double energy;
} reaction;

typedef struct {
	float xs[MAX_DATA_FILE];
	float ts[MAX_DATA_FILE];
	int maxIx;
	int ix;
} TimesVals;

typedef struct {
  int     maxEntries;   /* max number of entries in leaves */
  pcount *mwd_tree;
} mwdStore;

/* Linux: use the timeval struct to store seconds and microseconds
   Windows: use LARGE_INTEGER to store the stamps of QueryPerformanceCounter (divide by frequency later) */
typedef struct {
#if defined(__GNUC__)
	timeval start_t;
	timeval end_t;
#elif defined(_MSC_VER)
	LARGE_INTEGER start_t;
	LARGE_INTEGER end_t;
#endif
} Timer;

/* Windows: stores values of QueryPerformanceFrequency */ 
#if defined(_MSC_VER)
  LARGE_INTEGER frequency;
#endif

  typedef struct {
	chainLen *mols;
	pcount maxMolecules;
} MoleculeList;

typedef struct {
  pcount		no_of_mols;

  ptime			time;
  ptime			maxTime;
  ptime			nextSynchTime;
  ptime			synchTime;
  Timer			wallTime;

  long long events;
  long long nextSynchEvents;
  long long synchEvents;

#ifdef CALCMOMENTSOFDIST
  pcount momentDist[3];
#endif
#ifdef CALCFREEVOLUME
  double		freeVolumeFraction;
#endif

  float			conversion;

  int			noMoreReactions; // boolean flag set to true when no more events are possible

  reaction      reactions[NO_OF_REACTIONS];

  mwdStore      mwds[NO_OF_MOLSPECS][MAX_ARMS];
  int			arms[NO_OF_MOLSPECS];			// array of number of arms
  pcount        ms_cnts[NO_OF_MOLSPECS];
  pcount		react_cnts[NO_OF_REACTIONS];

  double 		volume;

  double		basetemp;
  double		deltatemp;
  double		temp;

  pcount 		initialMonomerMolecules;
  pcount 		currentMonomerMolecules;
  pcount		localMonomerParticles;
  
  probability   *scan_scratch;
  probability 	ratesVec[RATES_VEC_SIZE];
  int 			ratesVecPos;
  pcount 		tracerInitial;

  TimesVals 	timeCalcData;
  TimesVals		sysScaleData;
  float			scaleFactor;

  react_prob reactProbTree[2*REACT_PROB_TREE_LEAVES-1]; // was prob

  MoleculeList *expMols;
  
} sysState;


// Stores the header of the blocks of data that are communicated.
typedef struct {
	int 			stateTooBig;		// Flag indicating that size of 
										//    state has become too large.
	ptime           time;				// Time of system
	double			deltatemp;			// Temperature deviation of the system
	int 			noMoreReactions;	// Counts number of nodes that 
										//    have no events possible.
	pcount 			globalAllMonomer;	// Counts number of monomer 
										//    particles that have been 
										//    consumed and that still exist 
										//    as monomer.
} StatePacket;


void print_state();

inline const char *name(int index);
inline const char *rname(int index);
void print_reaction(int reaction_index);