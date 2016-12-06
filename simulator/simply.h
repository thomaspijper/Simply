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
#if defined(__GNUC__)
#include <sys/time.h>
#endif

 // Define inline directive for MSVC and GCC
#if defined(_MSC_VER)
  #define INLINE __inline
#elif defined(__GNUC__)
  #define INLINE inline
#endif

#define MAX_DATA_FILE 1000
typedef double probability;
typedef double ptime;
typedef unsigned long long pcount; // stores a number of particles
typedef unsigned long long rcount; // stores a number of reaction events
#ifdef LONGCHAINSUPPORT
typedef unsigned chainLen;
#define CHAINLENLIMIT UINT_MAX
#else
typedef short chainLen;
#define CHAINLENLIMIT SHRT_MAX
#endif

typedef struct {
	double	rc;
	int		arg_ms1;
	int		arg_ms2;
	int		res_ms1;
	int		res_ms2;    
	double	energy;
} reaction;

typedef struct {
	int		maxEntries;   /* max number of entries in leaves */
	pcount	*mwd_tree;
} mwdStore;

/* Linux: use the timeval struct to store seconds and microseconds
   Windows: use LARGE_INTEGER to store the stamps of QueryPerformanceCounter (divide by frequency later) */
typedef struct {
#if defined(__GNUC__)
	struct timeval start_t;
	struct timeval end_t;
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
	chainLen	*mols;
	pcount		maxMolecules;
} MoleculeList;

typedef struct {
	pcount			no_of_mols;

	ptime			time;
	ptime			maxTime;
	unsigned long long	nextSynchTime;	// In milliseconds. We cannot use a float type, since there is no guarantee that 0 can be represented exactly.
	unsigned long long	synchTime;		// In milliseconds. We cannot use a float type, since there is no guarantee that 0 can be represented exactly.
	Timer			wallTime;

	rcount			events;
	rcount			nextSynchEvents;
	rcount			synchEvents;

#ifdef CALCMOMENTSOFDIST
	pcount			momentDist[3];
#endif
#ifdef CALCFREEVOLUME
	double			freeVolumeFraction;
#endif

	float			conversion;

	int				noMoreReactionsLocal; // Boolean flag set to true when no more events are possible
	int				noMoreReactions; // Counts number of nodes that have no events possible

	reaction		reactions[NO_OF_REACTIONS];

	mwdStore		mwds[NO_OF_MOLSPECS][MAX_ARMS];
	int				arms[NO_OF_MOLSPECS]; // Array of number of arms
	pcount			ms_cnts[NO_OF_MOLSPECS];
	rcount			react_cnts[NO_OF_REACTIONS];

	double			volume;

	double			basetemp;
	double			deltatemp;
	double			temp;

	pcount 			initialMonomerMolecules;
	pcount 			currentMonomerMolecules;
	pcount			localMonomerParticles;

	double			scaleFactor;

	probability		reactProbTree[2 * REACT_PROB_TREE_LEAVES - 1]; // was prob

	MoleculeList	*expMols;
  
} sysState;


// Stores the header of the blocks of data that are communicated.
typedef struct {
	int 		stateTooBig;					// Flag indicating that size of state has become too large
	ptime		time;							// Time of system
	double		deltatemp;						// Temperature deviation of the system
	int 		noMoreReactions;				// Counts the number of nodes that have no events possible
	pcount 		globalAllMonomer;				// Counts the number of monomer particles that has been consumed + the number that still exists as monomer
	rcount		react_cnts[NO_OF_REACTIONS];	// Counts the number of reaction events per reaction
} StatePacket;


void print_state(void);
INLINE const char *name(int index);
INLINE const char *rname(int index);
void print_reaction(int reaction_index);
void print_kinetic_model(void);
void monomerAudit(const char *str);