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

#define GLOBAL_MONOMER_PARTICLES 496115969
#define SEED 0
#define CHANGESEED 0
#define MAX_SIM_TIME 300
#define MAX_WALL_TIME 30
#define MAX_EVENTS 0
#define MAX_CONVERSION 0.80F
#define SYNCH_TIME_INTERVAL 0
#define SYNCH_EVENTS_INTERVAL 5000000
#define MAXMONOMER 1
#define MAXSIMPLE 5
#define MAXPOLY 8
#define NO_OF_MOLSPECS 8
#define ARMS_INIT {1, 1, 1, 1, 1, 1, 1, 1}
#define MAX_ARMS 1
#define TOTAL_ARMS 8
#define M 0
#define RX 1
#define R_star 2
#define C 3
#define CX 4
#define D 5
#define P 6
#define PX 7
#define NO_MOL 8
#define POLY 9
#define SIMPLE 10
#define COMPLEX 11
#define init_activation 0
#define init_deactivation 1
#define initiation 2
#define propagation 3
#define activation 4
#define deactivation 5
#define disproportionation 6
#define combination 7
#define NO_OF_REACTIONS 8
#define MONOMERCONCENTRATION 9.350000000000000e+0
#define BASETEMP 343.0
#define STARTTEMP 343.0
#define SIMULATEHEATING 0
#define COOLINGRATE 0
#define RECALCCONVERSION 1
#define CALCMOMENTSOFDIST 1
#define CALCFREEVOLUME 1
#define SZYMANSKI 2
#define VF0 0.025
#define ALPHA_M 0.001
#define ALPHA_P 0.00048
#define TG_M 167.0
#define TG_P 378.0
#define LONGCHAINSUPPORT
#define MONO_AUDIT 1
#define EXPLICITSYSTEM 1
#define TREE_UPDATE_BODY {case init_activation:\
                            state.reactProbTree[7] = state.reactions[0].rc * (state.ms_cnts[RX] * state.ms_cnts[C]);\
                            state.reactProbTree[8] = state.reactions[1].rc * (state.ms_cnts[R_star] * state.ms_cnts[CX]);\
                            state.reactProbTree[9] = state.reactions[2].rc * (state.ms_cnts[R_star] * state.ms_cnts[M]);\
                            state.reactProbTree[11] = state.reactions[4].rc * (state.ms_cnts[PX] * state.ms_cnts[C]);\
                            state.reactProbTree[12] = state.reactions[5].rc * (state.ms_cnts[P] * state.ms_cnts[CX]);\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case init_deactivation:\
                            state.reactProbTree[7] = state.reactions[0].rc * (state.ms_cnts[RX] * state.ms_cnts[C]);\
                            state.reactProbTree[8] = state.reactions[1].rc * (state.ms_cnts[R_star] * state.ms_cnts[CX]);\
                            state.reactProbTree[9] = state.reactions[2].rc * (state.ms_cnts[R_star] * state.ms_cnts[M]);\
                            state.reactProbTree[11] = state.reactions[4].rc * (state.ms_cnts[PX] * state.ms_cnts[C]);\
                            state.reactProbTree[12] = state.reactions[5].rc * (state.ms_cnts[P] * state.ms_cnts[CX]);\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case initiation:\
                            state.reactProbTree[8] = state.reactions[1].rc * (state.ms_cnts[R_star] * state.ms_cnts[CX]);\
                            state.reactProbTree[9] = state.reactions[2].rc * (state.ms_cnts[R_star] * state.ms_cnts[M]);\
                            state.reactProbTree[10] = state.reactions[3].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[12] = state.reactions[5].rc * (state.ms_cnts[P] * state.ms_cnts[CX]);\
                            state.reactProbTree[13] = state.reactions[6].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[14] = state.reactions[7].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[6] = state.reactProbTree[13] + state.reactProbTree[14];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case propagation:\
                            state.reactProbTree[9] = state.reactions[2].rc * (state.ms_cnts[R_star] * state.ms_cnts[M]);\
                            state.reactProbTree[10] = state.reactions[3].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case activation:\
                            state.reactProbTree[7] = state.reactions[0].rc * (state.ms_cnts[RX] * state.ms_cnts[C]);\
                            state.reactProbTree[8] = state.reactions[1].rc * (state.ms_cnts[R_star] * state.ms_cnts[CX]);\
                            state.reactProbTree[10] = state.reactions[3].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[11] = state.reactions[4].rc * (state.ms_cnts[PX] * state.ms_cnts[C]);\
                            state.reactProbTree[12] = state.reactions[5].rc * (state.ms_cnts[P] * state.ms_cnts[CX]);\
                            state.reactProbTree[13] = state.reactions[6].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[14] = state.reactions[7].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[6] = state.reactProbTree[13] + state.reactProbTree[14];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case deactivation:\
                            state.reactProbTree[7] = state.reactions[0].rc * (state.ms_cnts[RX] * state.ms_cnts[C]);\
                            state.reactProbTree[8] = state.reactions[1].rc * (state.ms_cnts[R_star] * state.ms_cnts[CX]);\
                            state.reactProbTree[10] = state.reactions[3].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[11] = state.reactions[4].rc * (state.ms_cnts[PX] * state.ms_cnts[C]);\
                            state.reactProbTree[12] = state.reactions[5].rc * (state.ms_cnts[P] * state.ms_cnts[CX]);\
                            state.reactProbTree[13] = state.reactions[6].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[14] = state.reactions[7].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[6] = state.reactProbTree[13] + state.reactProbTree[14];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case disproportionation:\
                            state.reactProbTree[10] = state.reactions[3].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[12] = state.reactions[5].rc * (state.ms_cnts[P] * state.ms_cnts[CX]);\
                            state.reactProbTree[13] = state.reactions[6].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[14] = state.reactions[7].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[6] = state.reactProbTree[13] + state.reactProbTree[14];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;\
                          case combination:\
                            state.reactProbTree[10] = state.reactions[3].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                            state.reactProbTree[12] = state.reactions[5].rc * (state.ms_cnts[P] * state.ms_cnts[CX]);\
                            state.reactProbTree[13] = state.reactions[6].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[14] = state.reactions[7].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                            state.reactProbTree[6] = state.reactProbTree[13] + state.reactProbTree[14];\
                            state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                            state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                            state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                            state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                            state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];\
                            break;}
#define REACTION_PROBABILITY_TREE_INIT {state.reactProbTree[7] = state.reactions[0].rc * (state.ms_cnts[RX] * state.ms_cnts[C]);\
                                        state.reactProbTree[8] = state.reactions[1].rc * (state.ms_cnts[R_star] * state.ms_cnts[CX]);\
                                        state.reactProbTree[9] = state.reactions[2].rc * (state.ms_cnts[R_star] * state.ms_cnts[M]);\
                                        state.reactProbTree[10] = state.reactions[3].rc * (state.ms_cnts[P] * state.ms_cnts[M]);\
                                        state.reactProbTree[11] = state.reactions[4].rc * (state.ms_cnts[PX] * state.ms_cnts[C]);\
                                        state.reactProbTree[12] = state.reactions[5].rc * (state.ms_cnts[P] * state.ms_cnts[CX]);\
                                        state.reactProbTree[13] = state.reactions[6].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                                        state.reactProbTree[14] = state.reactions[7].rc * (state.ms_cnts[P] * (state.ms_cnts[P] - 1));\
                                        state.reactProbTree[6] = state.reactProbTree[13] + state.reactProbTree[14];\
                                        state.reactProbTree[5] = state.reactProbTree[11] + state.reactProbTree[12];\
                                        state.reactProbTree[4] = state.reactProbTree[9] + state.reactProbTree[10];\
                                        state.reactProbTree[3] = state.reactProbTree[7] + state.reactProbTree[8];\
                                        state.reactProbTree[2] = state.reactProbTree[5] + state.reactProbTree[6];\
                                        state.reactProbTree[1] = state.reactProbTree[3] + state.reactProbTree[4];\
                                        state.reactProbTree[0] = state.reactProbTree[1] + state.reactProbTree[2];}
#define REACT_PROB_TREE_LEAVES 8
#define REACTIONS_INIT {state.reactions[0].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (6.97316800e-08*1.88464000e+02*1) / ( (6.97316800e-08*1*exp(1.0/state.freeVolumeFraction)) + (1.88464000e+02*1));\
                        state.reactions[0].arg_ms1 = RX;\
                        state.reactions[0].arg_ms2 = C;\
                        state.reactions[0].res_ms1 = R_star;\
                        state.reactions[0].res_ms2 = CX;\
                        state.reactions[0].energy = 0;\
                        state.reactions[1].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1 / ((1/(9.95089920e-3)) + (1/((1.88464000e+02*1*exp(-1.0/state.freeVolumeFraction))+(1.97887200e-05*(state.ms_cnts[M]/AVOGADRO/state.volume))))));\
                        state.reactions[1].arg_ms1 = R_star;\
                        state.reactions[1].arg_ms2 = CX;\
                        state.reactions[1].res_ms1 = RX;\
                        state.reactions[1].res_ms2 = C;\
                        state.reactions[1].energy = 0;\
                        state.reactions[2].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1.98075664e-05*1.88464000e+08*1) / ( (1.98075664e-05*1*exp(1.0/state.freeVolumeFraction)) + (1.88464000e+08*1));\
                        state.reactions[2].arg_ms1 = R_star;\
                        state.reactions[2].arg_ms2 = M;\
                        state.reactions[2].res_ms1 = P;\
                        state.reactions[2].res_ms2 = NO_MOL;\
                        state.reactions[2].energy = 0;\
                        state.reactions[3].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1.98075664e-05*1.88464000e+08*1) / ( (1.98075664e-05*1*exp(1.0/state.freeVolumeFraction)) + (1.88464000e+08*1));\
                        state.reactions[3].arg_ms1 = P;\
                        state.reactions[3].arg_ms2 = M;\
                        state.reactions[3].res_ms1 = P;\
                        state.reactions[3].res_ms2 = NO_MOL;\
                        state.reactions[3].energy = 0;\
                        state.reactions[4].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (6.97316800e-08*1.88464000e+02*1) / ( (6.97316800e-08*1*exp(1.0/state.freeVolumeFraction)) + (1.88464000e+02*1));\
                        state.reactions[4].arg_ms1 = PX;\
                        state.reactions[4].arg_ms2 = C;\
                        state.reactions[4].res_ms1 = P;\
                        state.reactions[4].res_ms2 = CX;\
                        state.reactions[4].energy = 0;\
                        state.reactions[5].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1 / ((1/(9.95089920e-3)) + (1/((1.88464000e+02*1*exp(-1.0/state.freeVolumeFraction))+(1.97887200e-05*(state.ms_cnts[M]/AVOGADRO/state.volume))))));\
                        state.reactions[5].arg_ms1 = P;\
                        state.reactions[5].arg_ms2 = CX;\
                        state.reactions[5].res_ms1 = PX;\
                        state.reactions[5].res_ms2 = C;\
                        state.reactions[5].energy = 0;\
                        state.reactions[6].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1 / ((1/(3.76928000e-1)) + (1/((1.69994528e+06*((double)state.momentDist[0]/(double)state.momentDist[1])*((double)state.momentDist[0]/(double)state.momentDist[1])*exp(-1.0/state.freeVolumeFraction))+(6.79978112e+06*(state.ms_cnts[M]/AVOGADRO/state.volume))))));\
                        state.reactions[6].arg_ms1 = P;\
                        state.reactions[6].arg_ms2 = P;\
                        state.reactions[6].res_ms1 = D;\
                        state.reactions[6].res_ms2 = D;\
                        state.reactions[6].energy = 0;\
                        state.reactions[7].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1 / ((1/(3.47150688e+0)) + (1/((1.84694720e+05*((double)state.momentDist[0]/(double)state.momentDist[1])*((double)state.momentDist[0]/(double)state.momentDist[1])*exp(-1.0/state.freeVolumeFraction))+(7.38778880e+05*(state.ms_cnts[M]/AVOGADRO/state.volume))))));\
                        state.reactions[7].arg_ms1 = P;\
                        state.reactions[7].arg_ms2 = P;\
                        state.reactions[7].res_ms1 = D;\
                        state.reactions[7].res_ms2 = NO_MOL;\
                        state.reactions[7].energy = 0;}
#define RATES_UPDATE_BODY {state.reactions[0].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (6.97316800e-08*1.88464000e+02*1) / ( (6.97316800e-08*1*exp(1.0/state.freeVolumeFraction)) + (1.88464000e+02*1));\
                           state.reactions[1].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1 / ((1/(9.95089920e-3)) + (1/((1.88464000e+02*1*exp(-1.0/state.freeVolumeFraction))+(1.97887200e-05*(state.ms_cnts[M]/AVOGADRO/state.volume))))));\
                           state.reactions[2].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1.98075664e-05*1.88464000e+08*1) / ( (1.98075664e-05*1*exp(1.0/state.freeVolumeFraction)) + (1.88464000e+08*1));\
                           state.reactions[3].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1.98075664e-05*1.88464000e+08*1) / ( (1.98075664e-05*1*exp(1.0/state.freeVolumeFraction)) + (1.88464000e+08*1));\
                           state.reactions[4].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (6.97316800e-08*1.88464000e+02*1) / ( (6.97316800e-08*1*exp(1.0/state.freeVolumeFraction)) + (1.88464000e+02*1));\
                           state.reactions[5].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1 / ((1/(9.95089920e-3)) + (1/((1.88464000e+02*1*exp(-1.0/state.freeVolumeFraction))+(1.97887200e-05*(state.ms_cnts[M]/AVOGADRO/state.volume))))));\
                           state.reactions[6].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1 / ((1/(3.76928000e-1)) + (1/((1.69994528e+06*((double)state.momentDist[0]/(double)state.momentDist[1])*((double)state.momentDist[0]/(double)state.momentDist[1])*exp(-1.0/state.freeVolumeFraction))+(6.79978112e+06*(state.ms_cnts[M]/AVOGADRO/state.volume))))));\
                           state.reactions[7].rc = ((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * (1 / ((1/(3.47150688e+0)) + (1/((1.84694720e+05*((double)state.momentDist[0]/(double)state.momentDist[1])*((double)state.momentDist[0]/(double)state.momentDist[1])*exp(-1.0/state.freeVolumeFraction))+(7.38778880e+05*(state.ms_cnts[M]/AVOGADRO/state.volume))))));}
#define MOLECULENAMES {case M:\
                         return ("M");\
                       case RX:\
                         return ("RX");\
                       case R_star:\
                         return ("R_star");\
                       case C:\
                         return ("C");\
                       case CX:\
                         return ("CX");\
                       case D:\
                         return ("D");\
                       case P:\
                         return ("P");\
                       case PX:\
                         return ("PX");}
#define REACTIONNAMES {case init_activation:\
                         return ("init_activation");\
                       case init_deactivation:\
                         return ("init_deactivation");\
                       case initiation:\
                         return ("initiation");\
                       case propagation:\
                         return ("propagation");\
                       case activation:\
                         return ("activation");\
                       case deactivation:\
                         return ("deactivation");\
                       case disproportionation:\
                         return ("disproportionation");\
                       case combination:\
                         return ("combination");}
#define DO_REACT_BODY {case init_activation:\
                         no_of_res = 2;\
                         prod1_ind = R_star;\
                         prod2_ind = CX;\
                         break;\
                       case init_deactivation:\
                         no_of_res = 2;\
                         prod1_ind = RX;\
                         prod2_ind = C;\
                         break;\
                       case initiation:\
                         prod1_ind = P;\
                         prod1_lens[0] = 1;\
                         break;\
                       case propagation:\
                         prod1_ind = P;\
                         prod1_lens[0] = react1_lens[0] + 1;\
                         break;\
                       case activation:\
                         no_of_res = 2;\
                         prod1_ind = P;\
                         prod2_ind = CX;\
                         prod1_lens[0] = react1_lens[0];\
                         break;\
                       case deactivation:\
                         no_of_res = 2;\
                         prod1_ind = PX;\
                         prod2_ind = C;\
                         prod1_lens[0] = react1_lens[0];\
                         break;\
                       case disproportionation:\
                         no_of_res = 2;\
                         prod1_ind = D;\
                         prod2_ind = D;\
                         prod1_lens[0] = react1_lens[0];\
                         prod2_lens[0] = react2_lens[0];\
                         break;\
                       case combination:\
                         prod1_ind = D;\
                         prod1_lens[0] = react1_lens[0] + react2_lens[0];\
                         break;}
#define MWD_INITS {state.ms_cnts[M] = 496115969;\
                   state.ms_cnts[RX] = 2589354;\
                   state.ms_cnts[C] = 1294677;}
