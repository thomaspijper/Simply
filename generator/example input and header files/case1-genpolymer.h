/* 
 *
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
 *
 */

#include <math.h>

#define GLOBAL_MONOMER_PARTICLES 9992417015
#define CHANGESEED 0
#define MAX_SIM_TIME 300
#define MAX_WALL_TIME 30
#define MAX_EVENTS 0
#define MAX_CONVERSION 0.0F
#define SYNCH_TIME_INTERVAL 5000
#define SYNCH_EVENTS_INTERVAL 0
#define MAXMONOMER 1
#define MAXSIMPLE 6
#define MAXPOLY 15
#define NO_OF_MOLSPECS 18
#define ARMS_INIT {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2}
#define MAX_ARMS 2
#define TOTAL_ARMS 21
#define M 0
#define I 1
#define Junk 2
#define I_Star 3
#define RAFTR 4
#define R 5
#define QpreChain 6
#define QpreArm 7
#define A_coupled 8
#define Adead 9
#define D 10
#define P_RAFT 11
#define P 12
#define A_RAFT 13
#define A 14
#define Q 15
#define Qstar 16
#define QstarStar 17
#define NO_MOL 18
#define POLY 19
#define SIMPLE 20
#define COMPLEX 21
#define decomposition 0
#define decompositionNot 1
#define initiation 2
#define reInitiation 3
#define propagationP 4
#define propagationA 5
#define actAddP 6
#define actAddA 7
#define actFragPLeft 8
#define actFragPRight 9
#define actFragALeft 10
#define actFragARight 11
#define terminationAA 12
#define terminationAP 13
#define terminationPP 14
#define addRAFTPP 15
#define fragRAFT1PP 16
#define fragRAFT4PP 17
#define addRAFTAP 18
#define addRAFTAP1 19
#define fragRAFT1AP 20
#define fragRAFT2AP 21
#define addRAFTAA 22
#define fragRAFT1AA 23
#define fragRAFT2AA 24
#define NO_OF_REACTIONS 25
#define MONOMERCONCENTRATION 8.170000000000000e+0
#define BASETEMP 0.0
#define STARTTEMP 0.0
#define SIMULATEHEATING 0
#define COOLINGRATE 0
#define RECALCCONVERSION 0
#define CALCMOMENTSOFDIST 1
#define CALCFREEVOLUME 0
#define SZYMANSKI 2
#define VF0 -1
#define ALPHA_M -1
#define ALPHA_P -1
#define TG_M -1
#define TG_P -1
#define MONO_AUDIT 1
#define EXPLICITSYSTEM 1
#define TREE_UPDATE_BODY {case decomposition:\
                            state.reactProbTree[31] = state.reactions[0].rc * (state.ms_cnts[I]);\
                            state.reactProbTree[32] = state.reactions[1].rc * (state.ms_cnts[I]);\
                            state.reactProbTree[33] = state.reactions[2].rc * (state.ms_cnts[I_Star] * state.ms_cnts[M]);\
                            state.reactProbTree[16] = state.reactProbTree[33] + state.reactProbTree[34];\
                            state.reactProbTree[15] = state.reactProbTree[31] + state.reactProbTree[32];\
                            state.reactProbTree[7] = state.reactProbTree[15] + state.reactProbTree[16];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case decompositionNot:\
                            state.reactProbTree[31] = state.reactions[0].rc * (state.ms_cnts[I]);\
                            state.reactProbTree[32] = state.reactions[1].rc * (state.ms_cnts[I]);\
                            state.reactProbTree[15] = state.reactProbTree[31] + state.reactProbTree[32];\
                            state.reactProbTree[7] = state.reactProbTree[15] + state.reactProbTree[16];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case initiation:\
                            state.reactProbTree[33] = state.reactions[2].rc * (state.ms_cnts[I_Star] * state.ms_cnts[M]);\
                            state.reactProbTree[34] = state.reactions[3].rc * (state.ms_cnts[R] * state.ms_cnts[M]);\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[45] = state.reactions[14].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[16] = state.reactProbTree[33] + state.reactProbTree[34];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[7] = state.reactProbTree[15] + state.reactProbTree[16];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case reInitiation:\
                            state.reactProbTree[33] = state.reactions[2].rc * (state.ms_cnts[I_Star] * state.ms_cnts[M]);\
                            state.reactProbTree[34] = state.reactions[3].rc * (state.ms_cnts[R] * state.ms_cnts[M]);\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[43] = state.reactions[12].rc * (state.ms_cnts[A] * (state.ms_cnts[A] - 1));\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[16] = state.reactProbTree[33] + state.reactProbTree[34];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[7] = state.reactProbTree[15] + state.reactProbTree[16];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case propagationP:\
                            state.reactProbTree[33] = state.reactions[2].rc * (state.ms_cnts[I_Star] * state.ms_cnts[M]);\
                            state.reactProbTree[34] = state.reactions[3].rc * (state.ms_cnts[R] * state.ms_cnts[M]);\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[16] = state.reactProbTree[33] + state.reactProbTree[34];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[7] = state.reactProbTree[15] + state.reactProbTree[16];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case propagationA:\
                            state.reactProbTree[33] = state.reactions[2].rc * (state.ms_cnts[I_Star] * state.ms_cnts[M]);\
                            state.reactProbTree[34] = state.reactions[3].rc * (state.ms_cnts[R] * state.ms_cnts[M]);\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[16] = state.reactProbTree[33] + state.reactProbTree[34];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[7] = state.reactProbTree[15] + state.reactProbTree[16];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case actAddP:\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[39] = state.reactions[8].rc * (state.ms_cnts[QpreChain]);\
                            state.reactProbTree[40] = state.reactions[9].rc * (state.ms_cnts[QpreChain]);\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[45] = state.reactions[14].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[19] = state.reactProbTree[39] + state.reactProbTree[40];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[9] = state.reactProbTree[19] + state.reactProbTree[20];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case actAddA:\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[41] = state.reactions[10].rc * (state.ms_cnts[QpreArm]);\
                            state.reactProbTree[42] = state.reactions[11].rc * (state.ms_cnts[QpreArm]);\
                            state.reactProbTree[43] = state.reactions[12].rc * (state.ms_cnts[A] * (state.ms_cnts[A] - 1));\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[20] = state.reactProbTree[41] + state.reactProbTree[42];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[9] = state.reactProbTree[19] + state.reactProbTree[20];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case actFragPLeft:\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[39] = state.reactions[8].rc * (state.ms_cnts[QpreChain]);\
                            state.reactProbTree[40] = state.reactions[9].rc * (state.ms_cnts[QpreChain]);\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[45] = state.reactions[14].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[19] = state.reactProbTree[39] + state.reactProbTree[40];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[9] = state.reactProbTree[19] + state.reactProbTree[20];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case actFragPRight:\
                            state.reactProbTree[34] = state.reactions[3].rc * (state.ms_cnts[R] * state.ms_cnts[M]);\
                            state.reactProbTree[39] = state.reactions[8].rc * (state.ms_cnts[QpreChain]);\
                            state.reactProbTree[40] = state.reactions[9].rc * (state.ms_cnts[QpreChain]);\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[19] = state.reactProbTree[39] + state.reactProbTree[40];\
                            state.reactProbTree[16] = state.reactProbTree[33] + state.reactProbTree[34];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[9] = state.reactProbTree[19] + state.reactProbTree[20];\
                            state.reactProbTree[7] = state.reactProbTree[15] + state.reactProbTree[16];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case actFragALeft:\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[41] = state.reactions[10].rc * (state.ms_cnts[QpreArm]);\
                            state.reactProbTree[42] = state.reactions[11].rc * (state.ms_cnts[QpreArm]);\
                            state.reactProbTree[43] = state.reactions[12].rc * (state.ms_cnts[A] * (state.ms_cnts[A] - 1));\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[20] = state.reactProbTree[41] + state.reactProbTree[42];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[9] = state.reactProbTree[19] + state.reactProbTree[20];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case actFragARight:\
                            state.reactProbTree[34] = state.reactions[3].rc * (state.ms_cnts[R] * state.ms_cnts[M]);\
                            state.reactProbTree[41] = state.reactions[10].rc * (state.ms_cnts[QpreArm]);\
                            state.reactProbTree[42] = state.reactions[11].rc * (state.ms_cnts[QpreArm]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[20] = state.reactProbTree[41] + state.reactProbTree[42];\
                            state.reactProbTree[16] = state.reactProbTree[33] + state.reactProbTree[34];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[9] = state.reactProbTree[19] + state.reactProbTree[20];\
                            state.reactProbTree[7] = state.reactProbTree[15] + state.reactProbTree[16];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case terminationAA:\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[43] = state.reactions[12].rc * (state.ms_cnts[A] * (state.ms_cnts[A] - 1));\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case terminationAP:\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[43] = state.reactions[12].rc * (state.ms_cnts[A] * (state.ms_cnts[A] - 1));\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[45] = state.reactions[14].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case terminationPP:\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[45] = state.reactions[14].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case addRAFTPP:\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[45] = state.reactions[14].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[47] = state.reactions[16].rc * (state.ms_cnts[Q]);\
                            state.reactProbTree[48] = state.reactions[17].rc * (state.ms_cnts[Q]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[23] = state.reactProbTree[47] + state.reactProbTree[48];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case fragRAFT1PP:\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[45] = state.reactions[14].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[47] = state.reactions[16].rc * (state.ms_cnts[Q]);\
                            state.reactProbTree[48] = state.reactions[17].rc * (state.ms_cnts[Q]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[23] = state.reactProbTree[47] + state.reactProbTree[48];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case fragRAFT4PP:\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[45] = state.reactions[14].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[47] = state.reactions[16].rc * (state.ms_cnts[Q]);\
                            state.reactProbTree[48] = state.reactions[17].rc * (state.ms_cnts[Q]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[23] = state.reactProbTree[47] + state.reactProbTree[48];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case addRAFTAP:\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[43] = state.reactions[12].rc * (state.ms_cnts[A] * (state.ms_cnts[A] - 1));\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[51] = state.reactions[20].rc * (state.ms_cnts[Qstar]);\
                            state.reactProbTree[52] = state.reactions[21].rc * (state.ms_cnts[Qstar]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[25] = state.reactProbTree[51] + state.reactProbTree[52];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case addRAFTAP1:\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[45] = state.reactions[14].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[51] = state.reactions[20].rc * (state.ms_cnts[Qstar]);\
                            state.reactProbTree[52] = state.reactions[21].rc * (state.ms_cnts[Qstar]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[25] = state.reactProbTree[51] + state.reactProbTree[52];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case fragRAFT1AP:\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[43] = state.reactions[12].rc * (state.ms_cnts[A] * (state.ms_cnts[A] - 1));\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[51] = state.reactions[20].rc * (state.ms_cnts[Qstar]);\
                            state.reactProbTree[52] = state.reactions[21].rc * (state.ms_cnts[Qstar]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[25] = state.reactProbTree[51] + state.reactProbTree[52];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case fragRAFT2AP:\
                            state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[45] = state.reactions[14].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[51] = state.reactions[20].rc * (state.ms_cnts[Qstar]);\
                            state.reactProbTree[52] = state.reactions[21].rc * (state.ms_cnts[Qstar]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[25] = state.reactProbTree[51] + state.reactProbTree[52];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case addRAFTAA:\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[43] = state.reactions[12].rc * (state.ms_cnts[A] * (state.ms_cnts[A] - 1));\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[54] = state.reactions[23].rc * (state.ms_cnts[QstarStar]);\
                            state.reactProbTree[55] = state.reactions[24].rc * (state.ms_cnts[QstarStar]);\
                            state.reactProbTree[27] = state.reactProbTree[55] + state.reactProbTree[56];\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[13] = state.reactProbTree[27] + state.reactProbTree[28];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[6] = state.reactProbTree[13] + state.reactProbTree[14];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case fragRAFT1AA:\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[43] = state.reactions[12].rc * (state.ms_cnts[A] * (state.ms_cnts[A] - 1));\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[54] = state.reactions[23].rc * (state.ms_cnts[QstarStar]);\
                            state.reactProbTree[55] = state.reactions[24].rc * (state.ms_cnts[QstarStar]);\
                            state.reactProbTree[27] = state.reactProbTree[55] + state.reactProbTree[56];\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[13] = state.reactProbTree[27] + state.reactProbTree[28];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[6] = state.reactProbTree[13] + state.reactProbTree[14];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case fragRAFT2AA:\
                            state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                            state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                            state.reactProbTree[43] = state.reactions[12].rc * (state.ms_cnts[A] * (state.ms_cnts[A] - 1));\
                            state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                            state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                            state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                            state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                            state.reactProbTree[54] = state.reactions[23].rc * (state.ms_cnts[QstarStar]);\
                            state.reactProbTree[55] = state.reactions[24].rc * (state.ms_cnts[QstarStar]);\
                            state.reactProbTree[27] = state.reactProbTree[55] + state.reactProbTree[56];\
                            state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                            state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                            state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                            state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                            state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                            state.reactProbTree[13] = state.reactProbTree[27] + state.reactProbTree[28];\
                            state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                            state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                            state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                            state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                            state.reactProbTree[6] = state.reactProbTree[13] + state.reactProbTree[14];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;}
#define REACTION_PROBABILITY_TREE_INIT {state.reactProbTree[31] = state.reactions[0].rc * (state.ms_cnts[I]);\
                                        state.reactProbTree[32] = state.reactions[1].rc * (state.ms_cnts[I]);\
                                        state.reactProbTree[33] = state.reactions[2].rc * (state.ms_cnts[I_Star] * state.ms_cnts[M]);\
                                        state.reactProbTree[34] = state.reactions[3].rc * (state.ms_cnts[R] * state.ms_cnts[M]);\
                                        state.reactProbTree[35] = state.reactions[4].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                                        state.reactProbTree[36] = state.reactions[5].rc * (state.ms_cnts[A] * state.ms_cnts[M]);\
                                        state.reactProbTree[37] = state.reactions[6].rc * (state.ms_cnts[P] * state.ms_cnts[RAFTR]);\
                                        state.reactProbTree[38] = state.reactions[7].rc * (state.ms_cnts[A] * state.ms_cnts[RAFTR]);\
                                        state.reactProbTree[39] = state.reactions[8].rc * (state.ms_cnts[QpreChain]);\
                                        state.reactProbTree[40] = state.reactions[9].rc * (state.ms_cnts[QpreChain]);\
                                        state.reactProbTree[41] = state.reactions[10].rc * (state.ms_cnts[QpreArm]);\
                                        state.reactProbTree[42] = state.reactions[11].rc * (state.ms_cnts[QpreArm]);\
                                        state.reactProbTree[43] = state.reactions[12].rc * (state.ms_cnts[A] * (state.ms_cnts[A] - 1));\
                                        state.reactProbTree[44] = state.reactions[13].rc * (state.ms_cnts[A] * state.ms_cnts[P]);\
                                        state.reactProbTree[45] = state.reactions[14].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                                        state.reactProbTree[46] = state.reactions[15].rc * (state.ms_cnts[P] * state.ms_cnts[P_RAFT]);\
                                        state.reactProbTree[47] = state.reactions[16].rc * (state.ms_cnts[Q]);\
                                        state.reactProbTree[48] = state.reactions[17].rc * (state.ms_cnts[Q]);\
                                        state.reactProbTree[49] = state.reactions[18].rc * (state.ms_cnts[A] * state.ms_cnts[P_RAFT]);\
                                        state.reactProbTree[50] = state.reactions[19].rc * (state.ms_cnts[A_RAFT] * state.ms_cnts[P]);\
                                        state.reactProbTree[51] = state.reactions[20].rc * (state.ms_cnts[Qstar]);\
                                        state.reactProbTree[52] = state.reactions[21].rc * (state.ms_cnts[Qstar]);\
                                        state.reactProbTree[53] = state.reactions[22].rc * (state.ms_cnts[A] * state.ms_cnts[A_RAFT]);\
                                        state.reactProbTree[54] = state.reactions[23].rc * (state.ms_cnts[QstarStar]);\
                                        state.reactProbTree[55] = state.reactions[24].rc * (state.ms_cnts[QstarStar]);\
                                        state.reactProbTree[30] = state.reactProbTree[61] + state.reactProbTree[62];\
                                        state.reactProbTree[29] = state.reactProbTree[59] + state.reactProbTree[60];\
                                        state.reactProbTree[28] = state.reactProbTree[57] + state.reactProbTree[58];\
                                        state.reactProbTree[27] = state.reactProbTree[55] + state.reactProbTree[56];\
                                        state.reactProbTree[26] = state.reactProbTree[53] + state.reactProbTree[54];\
                                        state.reactProbTree[25] = state.reactProbTree[51] + state.reactProbTree[52];\
                                        state.reactProbTree[24] = state.reactProbTree[49] + state.reactProbTree[50];\
                                        state.reactProbTree[23] = state.reactProbTree[47] + state.reactProbTree[48];\
                                        state.reactProbTree[22] = state.reactProbTree[45] + state.reactProbTree[46];\
                                        state.reactProbTree[21] = state.reactProbTree[43] + state.reactProbTree[44];\
                                        state.reactProbTree[20] = state.reactProbTree[41] + state.reactProbTree[42];\
                                        state.reactProbTree[19] = state.reactProbTree[39] + state.reactProbTree[40];\
                                        state.reactProbTree[18] = state.reactProbTree[37] + state.reactProbTree[38];\
                                        state.reactProbTree[17] = state.reactProbTree[35] + state.reactProbTree[36];\
                                        state.reactProbTree[16] = state.reactProbTree[33] + state.reactProbTree[34];\
                                        state.reactProbTree[15] = state.reactProbTree[31] + state.reactProbTree[32];\
                                        state.reactProbTree[14] = state.reactProbTree[29] + state.reactProbTree[30];\
                                        state.reactProbTree[13] = state.reactProbTree[27] + state.reactProbTree[28];\
                                        state.reactProbTree[12] = state.reactProbTree[25] + state.reactProbTree[26];\
                                        state.reactProbTree[11] = state.reactProbTree[23] + state.reactProbTree[24];\
                                        state.reactProbTree[10] = state.reactProbTree[21] + state.reactProbTree[22];\
                                        state.reactProbTree[9] = state.reactProbTree[19] + state.reactProbTree[20];\
                                        state.reactProbTree[8] = state.reactProbTree[17] + state.reactProbTree[18];\
                                        state.reactProbTree[7] = state.reactProbTree[15] + state.reactProbTree[16];\
                                        state.reactProbTree[6] = state.reactProbTree[13] + state.reactProbTree[14];\
                                        state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                                        state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                                        state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                                        state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                                        state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                                        state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];}
#define REACT_PROB_TREE_LEAVES 32
#define REACTIONS_INIT {state.reactions[0].rc = 1.26080000e-04;\
                        state.reactions[0].arg_ms1 = I;\
                        state.reactions[0].arg_ms2 = NO_MOL;\
                        state.reactions[0].res_ms1 = I_Star;\
                        state.reactions[0].res_ms2 = I_Star;\
                        state.reactions[0].energy = 0;\
                        state.reactions[1].rc = 7.09200000e-05;\
                        state.reactions[1].arg_ms1 = I;\
                        state.reactions[1].arg_ms2 = NO_MOL;\
                        state.reactions[1].res_ms1 = Junk;\
                        state.reactions[1].res_ms2 = Junk;\
                        state.reactions[1].energy = 0;\
                        state.reactions[2].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 5.42899680e-07;\
                        state.reactions[2].arg_ms1 = I_Star;\
                        state.reactions[2].arg_ms2 = M;\
                        state.reactions[2].res_ms1 = P;\
                        state.reactions[2].res_ms2 = NO_MOL;\
                        state.reactions[2].energy = 0;\
                        state.reactions[3].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 5.42899680e-07;\
                        state.reactions[3].arg_ms1 = R;\
                        state.reactions[3].arg_ms2 = M;\
                        state.reactions[3].res_ms1 = A;\
                        state.reactions[3].res_ms2 = NO_MOL;\
                        state.reactions[3].energy = 0;\
                        state.reactions[4].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 5.42899680e-07;\
                        state.reactions[4].arg_ms1 = P;\
                        state.reactions[4].arg_ms2 = M;\
                        state.reactions[4].res_ms1 = P;\
                        state.reactions[4].res_ms2 = NO_MOL;\
                        state.reactions[4].energy = 0;\
                        state.reactions[5].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 5.42899680e-07;\
                        state.reactions[5].arg_ms1 = A;\
                        state.reactions[5].arg_ms2 = M;\
                        state.reactions[5].res_ms1 = A;\
                        state.reactions[5].res_ms2 = NO_MOL;\
                        state.reactions[5].energy = 0;\
                        state.reactions[6].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 8.17620000e-05;\
                        state.reactions[6].arg_ms1 = P;\
                        state.reactions[6].arg_ms2 = RAFTR;\
                        state.reactions[6].res_ms1 = QpreChain;\
                        state.reactions[6].res_ms2 = NO_MOL;\
                        state.reactions[6].energy = 0;\
                        state.reactions[7].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 8.17620000e-05;\
                        state.reactions[7].arg_ms1 = A;\
                        state.reactions[7].arg_ms2 = RAFTR;\
                        state.reactions[7].res_ms1 = QpreArm;\
                        state.reactions[7].res_ms2 = NO_MOL;\
                        state.reactions[7].energy = 0;\
                        state.reactions[8].rc = 1.00000000e+05;\
                        state.reactions[8].arg_ms1 = QpreChain;\
                        state.reactions[8].arg_ms2 = NO_MOL;\
                        state.reactions[8].res_ms1 = P;\
                        state.reactions[8].res_ms2 = RAFTR;\
                        state.reactions[8].energy = 0;\
                        state.reactions[9].rc = 1.00000000e+05;\
                        state.reactions[9].arg_ms1 = QpreChain;\
                        state.reactions[9].arg_ms2 = NO_MOL;\
                        state.reactions[9].res_ms1 = P_RAFT;\
                        state.reactions[9].res_ms2 = R;\
                        state.reactions[9].energy = 0;\
                        state.reactions[10].rc = 1.00000000e+05;\
                        state.reactions[10].arg_ms1 = QpreArm;\
                        state.reactions[10].arg_ms2 = NO_MOL;\
                        state.reactions[10].res_ms1 = A;\
                        state.reactions[10].res_ms2 = RAFTR;\
                        state.reactions[10].energy = 0;\
                        state.reactions[11].rc = 1.00000000e+05;\
                        state.reactions[11].arg_ms1 = QpreArm;\
                        state.reactions[11].arg_ms2 = NO_MOL;\
                        state.reactions[11].res_ms1 = A_RAFT;\
                        state.reactions[11].res_ms2 = R;\
                        state.reactions[11].energy = 0;\
                        state.reactions[12].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 1.63524000e-02;\
                        state.reactions[12].arg_ms1 = A;\
                        state.reactions[12].arg_ms2 = A;\
                        state.reactions[12].res_ms1 = A_coupled;\
                        state.reactions[12].res_ms2 = NO_MOL;\
                        state.reactions[12].energy = 0;\
                        state.reactions[13].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 4.08810000e-02;\
                        state.reactions[13].arg_ms1 = A;\
                        state.reactions[13].arg_ms2 = P;\
                        state.reactions[13].res_ms1 = Adead;\
                        state.reactions[13].res_ms2 = NO_MOL;\
                        state.reactions[13].energy = 0;\
                        state.reactions[14].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 1.63524000e-01;\
                        state.reactions[14].arg_ms1 = P;\
                        state.reactions[14].arg_ms2 = P;\
                        state.reactions[14].res_ms1 = D;\
                        state.reactions[14].res_ms2 = NO_MOL;\
                        state.reactions[14].energy = 0;\
                        state.reactions[15].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 4.08810000e-04;\
                        state.reactions[15].arg_ms1 = P;\
                        state.reactions[15].arg_ms2 = P_RAFT;\
                        state.reactions[15].res_ms1 = Q;\
                        state.reactions[15].res_ms2 = NO_MOL;\
                        state.reactions[15].energy = 0;\
                        state.reactions[16].rc = 1.00000000e+05;\
                        state.reactions[16].arg_ms1 = Q;\
                        state.reactions[16].arg_ms2 = NO_MOL;\
                        state.reactions[16].res_ms1 = P;\
                        state.reactions[16].res_ms2 = P_RAFT;\
                        state.reactions[16].energy = 0;\
                        state.reactions[17].rc = 1.00000000e+05;\
                        state.reactions[17].arg_ms1 = Q;\
                        state.reactions[17].arg_ms2 = NO_MOL;\
                        state.reactions[17].res_ms1 = P_RAFT;\
                        state.reactions[17].res_ms2 = P;\
                        state.reactions[17].energy = 0;\
                        state.reactions[18].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 4.08810000e-04;\
                        state.reactions[18].arg_ms1 = A;\
                        state.reactions[18].arg_ms2 = P_RAFT;\
                        state.reactions[18].res_ms1 = Qstar;\
                        state.reactions[18].res_ms2 = NO_MOL;\
                        state.reactions[18].energy = 0;\
                        state.reactions[19].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 4.08810000e-04;\
                        state.reactions[19].arg_ms1 = A_RAFT;\
                        state.reactions[19].arg_ms2 = P;\
                        state.reactions[19].res_ms1 = Qstar;\
                        state.reactions[19].res_ms2 = NO_MOL;\
                        state.reactions[19].energy = 0;\
                        state.reactions[20].rc = 1.00000000e+05;\
                        state.reactions[20].arg_ms1 = Qstar;\
                        state.reactions[20].arg_ms2 = NO_MOL;\
                        state.reactions[20].res_ms1 = A;\
                        state.reactions[20].res_ms2 = P_RAFT;\
                        state.reactions[20].energy = 0;\
                        state.reactions[21].rc = 1.00000000e+05;\
                        state.reactions[21].arg_ms1 = Qstar;\
                        state.reactions[21].arg_ms2 = NO_MOL;\
                        state.reactions[21].res_ms1 = A_RAFT;\
                        state.reactions[21].res_ms2 = P;\
                        state.reactions[21].energy = 0;\
                        state.reactions[22].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * 4.08810000e-04;\
                        state.reactions[22].arg_ms1 = A;\
                        state.reactions[22].arg_ms2 = A_RAFT;\
                        state.reactions[22].res_ms1 = QstarStar;\
                        state.reactions[22].res_ms2 = NO_MOL;\
                        state.reactions[22].energy = 0;\
                        state.reactions[23].rc = 1.00000000e+05;\
                        state.reactions[23].arg_ms1 = QstarStar;\
                        state.reactions[23].arg_ms2 = NO_MOL;\
                        state.reactions[23].res_ms1 = A;\
                        state.reactions[23].res_ms2 = A_RAFT;\
                        state.reactions[23].energy = 0;\
                        state.reactions[24].rc = 1.00000000e+05;\
                        state.reactions[24].arg_ms1 = QstarStar;\
                        state.reactions[24].arg_ms2 = NO_MOL;\
                        state.reactions[24].res_ms1 = A_RAFT;\
                        state.reactions[24].res_ms2 = A;\
                        state.reactions[24].energy = 0;}
#define RATES_UPDATE_BODY {}
#define MOLECULENAMES {case M:\
                         return ("M");\
                       case I:\
                         return ("I");\
                       case Junk:\
                         return ("Junk");\
                       case I_Star:\
                         return ("I_Star");\
                       case RAFTR:\
                         return ("RAFTR");\
                       case R:\
                         return ("R");\
                       case QpreChain:\
                         return ("QpreChain");\
                       case QpreArm:\
                         return ("QpreArm");\
                       case A_coupled:\
                         return ("A_coupled");\
                       case Adead:\
                         return ("Adead");\
                       case D:\
                         return ("D");\
                       case P_RAFT:\
                         return ("P_RAFT");\
                       case P:\
                         return ("P");\
                       case A_RAFT:\
                         return ("A_RAFT");\
                       case A:\
                         return ("A");\
                       case Q:\
                         return ("Q");\
                       case Qstar:\
                         return ("Qstar");\
                       case QstarStar:\
                         return ("QstarStar");}
#define REACTIONNAMES {case decomposition:\
                         return ("decomposition");\
                       case decompositionNot:\
                         return ("decompositionNot");\
                       case initiation:\
                         return ("initiation");\
                       case reInitiation:\
                         return ("reInitiation");\
                       case propagationP:\
                         return ("propagationP");\
                       case propagationA:\
                         return ("propagationA");\
                       case actAddP:\
                         return ("actAddP");\
                       case actAddA:\
                         return ("actAddA");\
                       case actFragPLeft:\
                         return ("actFragPLeft");\
                       case actFragPRight:\
                         return ("actFragPRight");\
                       case actFragALeft:\
                         return ("actFragALeft");\
                       case actFragARight:\
                         return ("actFragARight");\
                       case terminationAA:\
                         return ("terminationAA");\
                       case terminationAP:\
                         return ("terminationAP");\
                       case terminationPP:\
                         return ("terminationPP");\
                       case addRAFTPP:\
                         return ("addRAFTPP");\
                       case fragRAFT1PP:\
                         return ("fragRAFT1PP");\
                       case fragRAFT4PP:\
                         return ("fragRAFT4PP");\
                       case addRAFTAP:\
                         return ("addRAFTAP");\
                       case addRAFTAP1:\
                         return ("addRAFTAP1");\
                       case fragRAFT1AP:\
                         return ("fragRAFT1AP");\
                       case fragRAFT2AP:\
                         return ("fragRAFT2AP");\
                       case addRAFTAA:\
                         return ("addRAFTAA");\
                       case fragRAFT1AA:\
                         return ("fragRAFT1AA");\
                       case fragRAFT2AA:\
                         return ("fragRAFT2AA");}
#define DO_REACT_BODY {case decomposition:\
                         no_of_res = 2;\
                         prod1_ind = I_Star;\
                         prod2_ind = I_Star;\
                         break;\
                       case decompositionNot:\
                         no_of_res = 2;\
                         prod1_ind = Junk;\
                         prod2_ind = Junk;\
                         break;\
                       case initiation:\
                         prod1_ind = P;\
                         prod1_lens[0] = 1;\
                         break;\
                       case reInitiation:\
                         prod1_ind = A;\
                         prod1_lens[0] = 1;\
                         break;\
                       case propagationP:\
                         prod1_ind = P;\
                         prod1_lens[0] = react1_lens[0] + 1;\
                         break;\
                       case propagationA:\
                         prod1_ind = A;\
                         prod1_lens[0] = react1_lens[0] + 1;\
                         break;\
                       case actAddP:\
                         prod1_ind = QpreChain;\
                         prod1_lens[0] = react1_lens[0];\
                         break;\
                       case actAddA:\
                         prod1_ind = QpreArm;\
                         prod1_lens[0] = react1_lens[0];\
                         break;\
                       case actFragPLeft:\
                         no_of_res = 2;\
                         prod1_ind = P;\
                         prod2_ind = RAFTR;\
                         prod1_lens[0] = react1_lens[0];\
                         break;\
                       case actFragPRight:\
                         no_of_res = 2;\
                         prod1_ind = P_RAFT;\
                         prod2_ind = R;\
                         prod1_lens[0] = react1_lens[0];\
                         break;\
                       case actFragALeft:\
                         no_of_res = 2;\
                         prod1_ind = A;\
                         prod2_ind = RAFTR;\
                         prod1_lens[0] = react1_lens[0];\
                         break;\
                       case actFragARight:\
                         no_of_res = 2;\
                         prod1_ind = A_RAFT;\
                         prod2_ind = R;\
                         prod1_lens[0] = react1_lens[0];\
                         break;\
                       case terminationAA:\
                         prod1_ind = A_coupled;\
                         prod1_lens[0] = react1_lens[0] + react2_lens[0];\
                         break;\
                       case terminationAP:\
                         prod1_ind = Adead;\
                         prod1_lens[0] = react1_lens[0] + react2_lens[0];\
                         break;\
                       case terminationPP:\
                         prod1_ind = D;\
                         prod1_lens[0] = react1_lens[0] + react2_lens[0];\
                         break;\
                       case addRAFTPP:\
                         prod1_ind = Q;\
                         prod1_arms = 2;\
                         prod1_lens[0] = react1_lens[0];\
                         prod1_lens[1] = react2_lens[0];\
                         break;\
                       case fragRAFT1PP:\
                         no_of_res = 2;\
                         prod1_ind = P;\
                         prod2_ind = P_RAFT;\
                         prod1_lens[0] = react1_lens[0];\
                         prod2_lens[0] = react1_lens[1];\
                         break;\
                       case fragRAFT4PP:\
                         no_of_res = 2;\
                         prod1_ind = P_RAFT;\
                         prod2_ind = P;\
                         prod1_lens[0] = react1_lens[0];\
                         prod2_lens[0] = react1_lens[1];\
                         break;\
                       case addRAFTAP:\
                         prod1_ind = Qstar;\
                         prod1_arms = 2;\
                         prod1_lens[0] = react1_lens[0];\
                         prod1_lens[1] = react2_lens[0];\
                         break;\
                       case addRAFTAP1:\
                         prod1_ind = Qstar;\
                         prod1_arms = 2;\
                         prod1_lens[0] = react1_lens[0];\
                         prod1_lens[1] = react2_lens[0];\
                         break;\
                       case fragRAFT1AP:\
                         no_of_res = 2;\
                         prod1_ind = A;\
                         prod2_ind = P_RAFT;\
                         prod1_lens[0] = react1_lens[0];\
                         prod2_lens[0] = react1_lens[1];\
                         break;\
                       case fragRAFT2AP:\
                         no_of_res = 2;\
                         prod1_ind = A_RAFT;\
                         prod2_ind = P;\
                         prod1_lens[0] = react1_lens[0];\
                         prod2_lens[0] = react1_lens[1];\
                         break;\
                       case addRAFTAA:\
                         prod1_ind = QstarStar;\
                         prod1_arms = 2;\
                         prod1_lens[0] = react1_lens[0];\
                         prod1_lens[1] = react2_lens[0];\
                         break;\
                       case fragRAFT1AA:\
                         no_of_res = 2;\
                         prod1_ind = A;\
                         prod2_ind = A_RAFT;\
                         prod1_lens[0] = react1_lens[0];\
                         prod2_lens[0] = react1_lens[1];\
                         break;\
                       case fragRAFT2AA:\
                         no_of_res = 2;\
                         prod1_ind = A_RAFT;\
                         prod2_ind = A;\
                         prod1_lens[0] = react1_lens[0];\
                         prod2_lens[0] = react1_lens[1];\
                         break;}
#define MWD_INITS {state.ms_cnts[M] = 9992417015;\
                   state.ms_cnts[I] = 1467674;\
                   state.ms_cnts[RAFTR] = 6115310;}
