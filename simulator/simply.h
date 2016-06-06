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

#define MAX_DATA_FILE 1000
typedef double probability;
typedef double react_prob;

typedef int    length;
typedef double ptime;
typedef double rateCoefficient;
#ifdef LONGCHAINSUPPORT
typedef unsigned int chainLen;
#else
typedef unsigned short chainLen;
#endif

// stores a number of particles
typedef unsigned long long pcount; // 64-bit version
//typedef unsigned long pcount;  // 32-bit version

typedef int molecule_spec;

typedef struct {
  molecule_spec mspec;
  length        len; /* only relevant in case of Poly */  
} molecule;

typedef struct {
  rateCoefficient rc;
  molecule_spec arg_ms1;
  molecule_spec arg_ms2;
  molecule_spec res_ms1;
  molecule_spec res_ms2;    
  double energy;
} reaction;

typedef struct {
	float xs[MAX_DATA_FILE];
	float ts[MAX_DATA_FILE];
	int maxIx;
	int ix;
} TimesVals;

typedef struct {
  /* all arrays have the same length: one entry per
   * type of molecule in the system
   */
	int64_t *molCnts;
} specmap;

typedef struct {
  int     maxEntries;   /* max number of entries in leaves */
  pcount *mwd_tree;
} mwdStore;

#define RATES_VEC_SIZE 10

typedef struct {
	int64_t tv_sec;
	int64_t tv_usec;
} timeval;

typedef struct {
	timeval start_t;
	timeval end_t;
} Timer;

typedef struct {
	chainLen *mols;
	unsigned long long int maxMolecules;
} MoleculeList;

typedef struct {
  pcount		no_of_mols;

  ptime			time;
  ptime			maxTime;
  ptime			nextSynchTime;
  ptime			synchTime;
  Timer			wallTime;

  unsigned long long events;
  unsigned long long nextSynchEvents;
  unsigned long long synchEvents;

#ifdef CALCMOMENTSOFDIST
  unsigned long long zerothMoment;
  unsigned long long firstMoment;
  unsigned long long secondMoment;
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


#define MOLTYPE(MSPEC) \
  (MSPEC<MAXSIMPLE?SIMPLE:(MSPEC<MAXPOLY?POLY:(MSPEC<NO_MOL?COMPLEX:NO_MOL)))

#define SPECTYPE(MOLS) \
  (((MOLS)<MAXSIMPLE)?SIMPLE:(((MOLS)<MAXPOLY)?POLY:((MOLS)<NO_MOL)?COMPLEX:NO_MOL))



#define initMol(MOL,MOLSPEC,LEN) {	\
  (MOL).mspec   = MOLSPEC;	        \
  (MOL).len     = (LEN);}	        \


#define isPolyMol(MOL) 	                \
  (MOLTYPE(MOL) == POLY)                \

#define isComplexMol(MOL) 	                \
  (MOLTYPE(MOL) == COMPLEX)                \


#define isSimpleMol(MOL)                \
    (MOLTYPE(MOL) == SIMPLE)            \


#define isPolySpec(MOL) 	        \
  (SPECTYPE(MOL) == POLY)               \

#define isSimpleSpec(MOL)               \
    (MOLTYPE(MOL) == SIMPLE)            \

#define isComplexSpec(MOL)				\
	(MOLTYPE(MOL) == COMPLEX)


#define reactToSpecInd1(REACT)          \
  REACT.arg_ms1       

#define reactToSpecInd2(REACT)          \
  REACT.arg_ms2       

 
void print_state();

inline const char * name(int index);
inline const char * rname(int index);
void print_reaction(int reaction_index);


