# 
#  Copyright (c) 2016 Thomas Pijper
#
#  This file is part of Simply.
#
#  Simply is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 3 of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
# 


# Enable print() for Python 2.x
from __future__ import print_function

import argparse
import os
import sys
import linecache

from readInput import readInput
read = readInput()

from writeOutput import writeOutput
write = writeOutput()

# The beginning of the generated C header file: license and math.h include directive
header = """/* 
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
\n"""

# Visual Studio compiler needs this
_MCS_VER = """#if defined(_MSC_VER)
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
#endif
\n"""

def getOut(outputfile, deletionFlag):
    if deletionFlag == True and os.path.isfile(outputfile):
        os.remove(outputfile)
    sys.exit()

def main():
    # Start with defining the output file
    outputfile = 'genpolymer.h'

    # Get arguments
    parser = argparse.ArgumentParser(description='This program converts an Simply input file into a C header file')
    parser.add_argument('-i', metavar='input file', help='specifies the input file', nargs='?', default='!no_filename_has_been_provided!')
    parser.add_argument('-o', metavar='output file', help='specifies the output file (default = "genpolymer.h")', nargs='?', default='!no_filename_has_been_provided!')
    parser.add_argument('-f', help='forces the output file to be overwritten', action='store_true')
    parser.add_argument('-v', help='prints detailed information', action='store_true')
    args = parser.parse_args()

    # Set the name of the output file
    if (args.o != '!no_filename_has_been_provided!'):
        outputfile = args.o

    # Set the name of the input file. Work around raw_input() -> input()
    version = sys.version_info
    if (args.i == '!no_filename_has_been_provided!'):
        if version[0] == 2:
            inputfile = raw_input('\nWhat is the name of your input file? ')
        elif version[0] == 3:
            inputfile = input('\nWhat is the name of your input file? ')
        else:
            print('\nCould not recognize the version of Python\n')
            sys.exit()
        if inputfile == '':
            print('\nNo filename has been provided. Exiting...')
            sys.exit()
    else:
        inputfile = args.i
    
    # Verify whether the input file exists
    if os.path.isfile(inputfile) == False:
        print('\nThe file you specified could not be found. Exiting...')
        sys.exit()

    # Retreive parameters while performing input checks
    generalDict = read.readGeneral(inputfile)
    if generalDict == -1:
        getOut(outputfile, args.f)
    moleculesList = read.readMolecules(inputfile, generalDict)
    if moleculesList == -1:
        getOut(outputfile, args.f)
    ratesList = read.readRateCoeffs(inputfile, generalDict, moleculesList)
    if ratesList == -1:
        getOut(outputfile, args.f)
    reactionsList = read.readReactions(inputfile, moleculesList, ratesList)
    if reactionsList == -1:
        getOut(outputfile, args.f)
    enthalpiesList = read.readEnthalpies(inputfile, generalDict, reactionsList)
    if enthalpiesList == -1:
        getOut(outputfile, args.f)
    freeVolumeDict = read.readFreeVolume(inputfile, generalDict, ratesList)
    if freeVolumeDict == -1:
        getOut(outputfile, args.f)

    if args.v == True:
        print('\nGeneral input:')
        print(generalDict)
        print('\nMolecules:')
        for i in range(len(moleculesList)):
            print(moleculesList[i])
        print('\nRate coefficients:')
        for i in range(len(ratesList)):
            print(ratesList[i])
        print('\nReactions:')
        for i in range(len(reactionsList)):
            print(reactionsList[i])
        print('\nEnthalpies:')
        for i in range(len(enthalpiesList)):
            print(enthalpiesList[i])
        print('\nFree volume parameters:')
        print(freeVolumeDict)

    # Check if output file exists
    if os.path.isfile(outputfile) == True and args.f == False:
        print('\nFile {0} already exists. Exiting...'.format(outputfile))
        sys.exit()
    elif os.path.isfile(outputfile) == True and args.f == True:
        os.remove(outputfile)
    print('\nWriting to file '+outputfile)

    # Write output
    write.header(outputfile,header)
    write._MCS_VER(outputfile,_MCS_VER)
    write.GLOBAL_MONOMER_PARTICLES(outputfile,generalDict,moleculesList)
    write.SEED(outputfile,generalDict)
    write.endingCriteria(outputfile,generalDict)
    write.syncMode(outputfile,generalDict)
    write.MAXMONOMER(outputfile,generalDict)    
    write.MAXSIMPLE(outputfile,moleculesList)
    write.MAXPOLY(outputfile,moleculesList)
    write.NO_OF_MOLSPECS(outputfile,moleculesList)
    write.arms(outputfile,moleculesList)
    write.molecules(outputfile,moleculesList)
    write.reactions(outputfile,reactionsList)
    write.NO_OF_REACTIONS(outputfile,reactionsList)
    write.MONOMERCONCENTRATION(outputfile,generalDict,moleculesList)
    write.BASETEMP(outputfile,generalDict)
    write.STARTTEMP(outputfile,generalDict)
    write.SIMULATEHEATING(outputfile,generalDict)
    write.COOLINGRATE(outputfile,generalDict)
    write.RECALCCONVERSION(outputfile,generalDict)
    write.CALCMOMENTSOFDIST(outputfile,generalDict,ratesList)
    write.CALCFREEVOLUME(outputfile,generalDict)
    write.SZYMANSKI(outputfile,generalDict)
    write.freeVolumeParameters(outputfile,generalDict,freeVolumeDict)
    write.LONGCHAINSUPPORT(outputfile,generalDict)
    write.TREE_UPDATE_BODY(outputfile,reactionsList,moleculesList)
    write.REACTION_PROBABILITY_TREE_INIT(outputfile,reactionsList,moleculesList)
    write.REACT_PROB_TREE_LEAVES(outputfile,reactionsList)
    write.REACTIONS_INIT(outputfile,generalDict,moleculesList, reactionsList, ratesList, enthalpiesList)
    write.RATES_UPDATE_BODY(outputfile,generalDict,moleculesList, reactionsList, ratesList)
    write.MOLECULENAMES(outputfile,moleculesList)
    write.REACTIONNAMES(outputfile,reactionsList)
    write.DO_REACT_BODY(outputfile,generalDict,moleculesList,reactionsList)
    write.MWD_INITS(outputfile,generalDict,moleculesList)

if __name__ == '__main__':
    # Execute the program
    main()
