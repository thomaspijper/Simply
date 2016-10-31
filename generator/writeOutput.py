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

import sys
from decimal import *

from readInput import readInput
read = readInput()

# For unit testing
try:
    from StringIO import StringIO  # For Python 2.x
except ImportError:
    from io import StringIO  # For Python 3.x

class writeOutput(object):

    def __init__(self):
        self.gasconstant = 8.3144598
        self.Avogadro = 6.022140857e23


    def header(self, file, header):
    # Write the header with license details
        self.writeSingleString(file, header)
        return True


    def _MCS_VER(self, file, _MCS_VER):
    # Write the header with license details
        self.writeSingleString(file, _MCS_VER)
        return True


    def GLOBAL_MONOMER_PARTICLES(self, file, generalDict, moleculesList):
    # Calculate the number of monomer particles and write it to the file

        particlecount = generalDict['particlecount']
        monomerList = generalDict['monomernames']
        concentration = self.calcConcentration(moleculesList)
        count = 0

        for i in range(len(monomerList)):
            for j in range(len(moleculesList)):
                if monomerList[i] == moleculesList[j][1]:
                    count = count + int(round(moleculesList[j][2] * particlecount / concentration))
                    break

        line = '#define GLOBAL_MONOMER_PARTICLES {0}\n'.format(count)
        self.writeSingleString(file, line)

        return True


    def SEED(self, file, generalDict):
    # Write the seed value (if specified)

        seedValue = generalDict['seed']

        if seedValue == -1:
            # Don't write anything, a seed will be automatically chosen
            return True

        line = '#define SEED {0}\n'.format(seedValue)
        self.writeSingleString(file, line)
        return True


    def endingCriteria(self, file, generalDict):
    # Write ending criteria

        line = '#define MAX_SIM_TIME {0}\n'.format(generalDict['maxsimtime'])
        self.writeSingleString(file, line)
        line = '#define MAX_WALL_TIME {0}\n'.format(generalDict['maxwalltime'])
        self.writeSingleString(file, line)
        line = '#define MAX_EVENTS {0}\n'.format(generalDict['maxevents'])
        self.writeSingleString(file, line)
        line = '#define MAX_CONVERSION {0}F\n'.format(round(generalDict['maxconversion'],5))
        self.writeSingleString(file, line)
        return True


    def syncMode(self, file, generalDict):
    # Write sync parameters

        line = '#define SYNCH_TIME_INTERVAL {0}\n'.format(generalDict['syncsimtime'])
        self.writeSingleString(file, line)
        line = '#define SYNCH_EVENTS_INTERVAL {0}\n'.format(generalDict['syncevents'])
        self.writeSingleString(file, line)
        return True


    def MAXMONOMER(self, file, generalDict):
    # Determine the number of different monomer species and write it to the file
        line = '#define MAXMONOMER {0}\n'.format(len(generalDict['monomernames']))
        self.writeSingleString(file, line)
        return True


    def MAXSIMPLE(self, file, moleculesList):
    # Determine the number of 'simple' molecular species and write it to the file
        counter = 0
        for i in range(len(moleculesList)):
            if moleculesList[i][0] == 'simple':
                counter = counter + 1

        line = '#define MAXSIMPLE {0}\n'.format(counter)
        self.writeSingleString(file, line)
        return True


    def MAXPOLY(self, file, moleculesList):
    # Determine the total number of 'simple' and 'poly' molecular species and write it to the file
        counter = 0
        for i in range(len(moleculesList)):
            if moleculesList[i][0] in ['simple','poly']:
                counter = counter + 1

        line = '#define MAXPOLY {0}\n'.format(counter)
        self.writeSingleString(file, line)
        return True


    def NO_OF_MOLSPECS(self, file, moleculesList):
    # Determine the total number of molecular species and write it to the file
        line = '#define NO_OF_MOLSPECS {0}\n'.format(len(moleculesList))
        self.writeSingleString(file, line)
        return True


    def ARMS_INIT(self, file, moleculesList):
    # Compile a list containing the maximum number of arms for each species
        list = []
        for i in range(len(moleculesList)):
            if moleculesList[i][0] in ['simple','poly']:
                list.append(1)
            elif moleculesList[i][0] == 'complex':
                list.append(moleculesList[i][2])

        line = '#define ARMS_INIT {{{0}}}\n'.format(str(list).strip('[]'))
        self.writeSingleString(file, line)
        return True


    def molecules(self, file, moleculesList):
    # Print a line for each molecular species, end one line for each general type
        number = 0
        for i in range(len(moleculesList)):
            line = '#define {0} {1}\n'.format(moleculesList[i][1], number)
            self.writeSingleString(file, line)
            number = number + 1
        line = '#define {0} {1}\n'.format('NO_MOL', number)
        self.writeSingleString(file, line)
        line = '#define {0} {1}\n'.format('POLY', number+1)
        self.writeSingleString(file, line)
        line = '#define {0} {1}\n'.format('SIMPLE', number+2)
        self.writeSingleString(file, line)
        line = '#define {0} {1}\n'.format('COMPLEX', number+3)
        self.writeSingleString(file, line)
        return True        


    def reactions(self, file, reactionsList):
    # Print a line for each reaction
        number = 0
        for i in range(len(reactionsList)):
            line = '#define {0} {1}\n'.format(reactionsList[i][0], number)
            self.writeSingleString(file, line)
            number = number + 1
        return True        


    def NO_OF_REACTIONS(self, file, reactionsList):
    # Determine the total number of reactions and write it to the file
        line = '#define NO_OF_REACTIONS {0}\n'.format(len(reactionsList))
        self.writeSingleString(file, line)
        return True


    def CONCENTRATION(self, file, generalDict, moleculesList):
    # Write the starting temperature of the reaction mixture to the file

        count = 0
        monomerList = generalDict['monomernames']
        for i in range(len(monomerList)):
            for j in range(len(moleculesList)):
                if monomerList[i] == moleculesList[j][1]:
                    count = count + moleculesList[j][2]
                    break

        line = '#define CONCENTRATION {0:.8e}\n'.format(Decimal(count))
        self.writeSingleString(file, line)
        return True


    def BASETEMP(self, file, generalDict):
    # Write the concentration of the reaction mixture to the file

        line = '#define BASETEMP {temperature}\n'.format(**generalDict)
        self.writeSingleString(file, line)
        return True


    def STARTTEMP(self, file, generalDict):
    # Write the concentration of the reaction mixture to the file
        
        if generalDict['coolingrate'] == 0.0:  # No cooling, so we set the starttemp equal to the basetemp
            line = '#define STARTTEMP {temperature}\n'.format(**generalDict)
        else:
            line = '#define STARTTEMP {starttemperature}\n'.format(**generalDict)
        self.writeSingleString(file, line)
        return True


    def SIMULATEHEATING(self, file, generalDict):
    # Write value indicating we want to simulate heating

        if generalDict['simulateheating'] != 0:
            line = '#define SIMULATEHEATING\n'
            self.writeSingleString(file, line)
        return True


    def COOLINGRATE(self, file, generalDict):
    # Write the cooling rate of the system to the file

        if (generalDict['coolingrate'] != 0):
            line = '#define COOLINGRATE {coolingrate}\n'.format(**generalDict)
            self.writeSingleString(file, line)
        return True


    def RECALCCONVERSION(self, file, generalDict):
    # Indicate whether the conversion should be recalculated on each iteration

        if (generalDict['recalcconversion'] != 0):
            line = '#define RECALCCONVERSION\n'
            self.writeSingleString(file, line)
        return True


    def CALCMOMENTSOFDIST(self, file, generalDict, ratesList):
    # Indicate whether the conversion should be recalculated on each iteration

        for i in range(len(ratesList)):
            if ratesList[i][1] in ['fp','ap','fpr','apr']:
                line = '#define CALCMOMENTSOFDIST\n'
                self.writeSingleString(file, line)
                break
            if i == (len(ratesList) - 1):
                if generalDict['calcdist'] != 0:
                    line = '#define CALCMOMENTSOFDIST\n'
                    self.writeSingleString(file, line)
        return True


    def CALCFREEVOLUME(self, file, generalDict):
    # Indicate whether the free volume fraction should be calculated

        if generalDict['freevolume'] != 0:
            line = '#define CALCFREEVOLUME\n'
            self.writeSingleString(file, line)
        return True


    def freeVolumeParameters(self, file, generalDict, freeVolumeDict):
    # Write the free volume parameters

        if generalDict['freevolume'] != 0:
            line = '#define VF0 {vf0}\n'.format(**freeVolumeDict)
            line = line + '#define ALPHA_M {alpham}\n'.format(**freeVolumeDict)
            line = line + '#define ALPHA_P {alphap}\n'.format(**freeVolumeDict)
            line = line + '#define TG_M {tgm}\n'.format(**freeVolumeDict)
            line = line + '#define TG_P {tgp}\n'.format(**freeVolumeDict)
            self.writeSingleString(file, line)
        return True

        
    def LONGCHAINSUPPORT(self, file, generalDict):
    # Indicate whether a 'short' or an 'int' should be used for storing chain lengths

        if generalDict['longchainsupport'] == 1:
            line = '#define LONGCHAINSUPPORT\n'
            self.writeSingleString(file, line)
        return True
        

    def REACTION_PROBABILITY_TREE_INIT(self, file, reactionsList, moleculesList):
    # We're going to define the reaction probability tree
    
        startFirstLine = '#define REACTION_PROBABILITY_TREE_INIT {'
        
        # First calculate the number of probability tree leaves
        i = 1
        while True:
            leavesNo = 2 ** i
            if leavesNo >= len(reactionsList):
                break
            i = i + 1

        # Write a line for each reaction leave
        for i in range(len(reactionsList)):

            # Write the leave index and reference to the reaction coefficient
            x = leavesNo - 1 + i
            if i == 0:
                line = startFirstLine + 'state.reactProbTree[{0}] = state.reactions[{1}].rc * '.format(x,i)
            else:
                line = len(startFirstLine) * ' ' + 'state.reactProbTree[{0}] = state.reactions[{1}].rc * '.format(x,i)
            self.writeSingleString(file, line)

            # Write the reactant(s)
            if reactionsList[i][3] == '':
                reactant1 = read.isReactingComplex('',reactionsList[i][2],moleculesList)[1]
                line = '(state.ms_cnts[{0}]);'.format(reactant1)
                self.writeSingleString(file, line)
            elif reactionsList[i][2] == reactionsList[i][3]:
                reactant1 = read.isReactingComplex('',reactionsList[i][2],moleculesList)[1]
                reactant2 = read.isReactingComplex('',reactionsList[i][3],moleculesList)[1]
                line = '(state.ms_cnts[{0}] * (state.ms_cnts[{0}] - 1));'.format(reactant1,reactant2)
                self.writeSingleString(file, line)
            else:
                reactant1 = read.isReactingComplex('',reactionsList[i][2],moleculesList)[1]
                reactant2 = read.isReactingComplex('',reactionsList[i][3],moleculesList)[1]
                line = '(state.ms_cnts[{0}] * state.ms_cnts[{1}]);'.format(reactant1,reactant2)
                self.writeSingleString(file, line)

            # Finish the line
            line = '\\\n'
            self.writeSingleString(file, line)

        # Write a line for each reaction branch
        for i in range(leavesNo - 2, -1, -1):
            line = len(startFirstLine) * ' ' + 'state.reactProbTree[{0}] = state.reactProbTree[{1}] + state.reactProbTree[{2}];'.format(i, i*2+1, i*2+2)
            self.writeSingleString(file, line)

            # Finish the line
            if i == 0:     # We're counting down
                line = '}\n'
            else:
                line = '\\\n'
            self.writeSingleString(file, line)

        return True


    def TREE_UPDATE_BODY(self, file, reactionsList, moleculesList):
    # We're going to define an update of the reaction probability tree for each reaction case
    
        # First calculate the number of probability tree leaves
        exp = 1
        while True:
            leavesNo = 2 ** exp
            if leavesNo >= len(reactionsList):
                break
            exp += 1

        # Write the beginning of the first line
        startFirstLine = '#define TREE_UPDATE_BODY {'
        self.writeSingleString(file, startFirstLine)

        # We're going to make a case for each reaction
        for i in range(len(reactionsList)):
            indicesAffectedReactions = [] # For storing a list of reactions that require updating
            
            # Write the first line
            if i == 0:
                line = 'case {0}:\\\n{1}'.format(reactionsList[i][0], len(startFirstLine)*' ')
                self.writeSingleString(file, line)
            else:
                line = len(startFirstLine) * ' ' + 'case {0}:\\\n{1}'.format(reactionsList[i][0], len(startFirstLine)*' ')
                self.writeSingleString(file, line)
            
            # Find out which molcounts are affected by this reaction
            # Compile list of reactants and products
            molChanged = [read.isReactingComplex('',reactionsList[i][2],moleculesList)[1],
            read.isReactingComplex('',reactionsList[i][3],moleculesList)[1],
            read.isReactingComplex('',reactionsList[i][4],moleculesList)[1],
            read.isReactingComplex('',reactionsList[i][5],moleculesList)[1]]
            # Discard items that do not appear or disappear (for (de)propagations)
            if molChanged[0] == molChanged[2]:
                del molChanged[2]
                del molChanged[0]
            molChanged[:] = [x for x in molChanged if x != ''] # Remove empty items

            # Cycle through the list of reactions and find which reactions are affected
            for j in range(len(reactionsList)):
                r1 = read.isReactingComplex('',reactionsList[j][2],moleculesList)[1]
                r2 = read.isReactingComplex('',reactionsList[j][3],moleculesList)[1]
                r = [r1,r2] # The reactants of this reaction
                if bool(set(molChanged) & set(r)) == True:
                    indicesAffectedReactions.append(leavesNo - 1 + j)
            
                    # Write the leave index and reaction coefficient
                    x = leavesNo - 1 + j
                    line = '  state.reactProbTree[{0}] = state.reactions[{1}].rc * '.format(x,j)
                    self.writeSingleString(file, line)
                    # Write the reactant(s)
                    if reactionsList[j][3] == '':
                        reactant1 = read.isReactingComplex('',reactionsList[j][2],moleculesList)[1]
                        line = '(state.ms_cnts[{0}]);\\\n'.format(reactant1)
                        self.writeSingleString(file, line)
                    elif reactionsList[j][2] == reactionsList[j][3]:
                        reactant1 = read.isReactingComplex('',reactionsList[j][2],moleculesList)[1]
                        reactant2 = read.isReactingComplex('',reactionsList[j][3],moleculesList)[1]
                        line = '(state.ms_cnts[{0}] * (state.ms_cnts[{0}] - 1));\\\n'.format(reactant1,reactant2)
                        self.writeSingleString(file, line)
                    else:
                        reactant1 = read.isReactingComplex('',reactionsList[j][2],moleculesList)[1]
                        reactant2 = read.isReactingComplex('',reactionsList[j][3],moleculesList)[1]
                        line = '(state.ms_cnts[{0}] * state.ms_cnts[{1}]);\\\n'.format(reactant1,reactant2)
                        self.writeSingleString(file, line)
                    # Write some whitespace
                    line = len(startFirstLine) * ' '
                    self.writeSingleString(file, line)

            # Now, we're going to write a line for each reaction branch
            indicesAffectedReactions = list(set(indicesAffectedReactions)) # Remove duplicates of indices
            indicesAffectedReactions.sort(reverse=True)

            indicesAffectedBranches = []            
            for j in range(len(indicesAffectedReactions)):
                x = indicesAffectedReactions[j]
                while (x != 0):
                    if x % 2 == 0:
                        x = x - 1
                    x = (x - 1) / 2
                    indicesAffectedBranches.append(int(x))

            indicesAffectedBranches = list(set(indicesAffectedBranches)) # Remove duplicates of indices
            indicesAffectedBranches.sort(reverse=True)

            for k in range(len(indicesAffectedBranches)):
                y = indicesAffectedBranches[k]
                line = '  state.reactProbTree[{0}] = state.reactProbTree[{1}] + state.reactProbTree[{2}];\\\n'.format(y, y*2+1, y*2+2)
                self.writeSingleString(file, line)
                
                # Write some whitespace
                line = len(startFirstLine)*' '
                self.writeSingleString(file, line)

            if i == (len(reactionsList) - 1):
                line = '  break;}\n'
                self.writeSingleString(file, line)
            else:
                line = '  break;\\\n'
                self.writeSingleString(file, line)

        return True          


    def REACT_PROB_TREE_LEAVES(self, file, reactionsList):
    # Write the number of leaves

        i = 1
        while True:
            leavesNo = 2 ** i
            if leavesNo >= len(reactionsList):
                break
            i = i + 1

        line = '#define REACT_PROB_TREE_LEAVES {0}\n'.format(leavesNo)
        self.writeSingleString(file, line)
        return True
 
   
    def MAX_ARMS(self, file, reactionsList):
    # Write the maximum number of arms

        maxArms = 1
        for i in range(len(reactionsList)):
            if reactionsList[i][0] == 'complex' and reactionsList[i][2] > maxArms:
                maxArms = reactionsList[i][2]

        line = '#define MAX_ARMS {0}\n'.format(maxArms)
        self.writeSingleString(file, line)
        return True
 

    def REACTIONS_INIT(self, file, generalDict, moleculesList, reactionsList, ratesList, enthalpiesList):
    # We're going to define reactants, products, and rate coefficient for each reaction

        particlecount = generalDict['particlecount']
        concentration = self.calcConcentration(moleculesList)
        startFirstLine = '#define REACTIONS_INIT {'

        # Go through the list of reactions
        for i in range(len(reactionsList)):

            # Look up the corresponding rate coefficient details
            for j in range(len(ratesList)):
                if reactionsList[i][1] == ratesList[j][0]:
                    rate = ratesList[j]
                    break

            # Compile the string that describes the rate coefficient
            rcline = self.compileRateString(rate,particlecount,concentration,reactionsList[i],i) + '\\\n'

            # Next, get the reactants and products
            r1 = read.isReactingComplex('', reactionsList[i][2], moleculesList)[1]
            if reactionsList[i][3] == '':
                r2 = 'NO_MOL'
            else:
                r2 = read.isReactingComplex('', reactionsList[i][3], moleculesList)[1]
            p1 = read.isReactingComplex('', reactionsList[i][4], moleculesList)[1]
            if reactionsList[i][5] == '':
                p2 = 'NO_MOL'
            else:
                p2 = read.isReactingComplex('', reactionsList[i][5], moleculesList)[1]

            # Write rate coefficients, reactants, and products to header file
            if i == 0:
                line = startFirstLine + rcline
            else:
                line = len(startFirstLine) * ' ' +  rcline
            self.writeSingleString(file, line)
            line = len(startFirstLine) * ' ' + 'state.reactions[{0}].arg_ms1 = {1};\\\n'.format(i, r1)
            self.writeSingleString(file, line)
            line = len(startFirstLine) * ' ' + 'state.reactions[{0}].arg_ms2 = {1};\\\n'.format(i, r2)
            self.writeSingleString(file, line)
            line = len(startFirstLine) * ' ' + 'state.reactions[{0}].res_ms1 = {1};\\\n'.format(i, p1)
            self.writeSingleString(file, line)
            line = len(startFirstLine) * ' ' + 'state.reactions[{0}].res_ms2 = {1};'.format(i, p2)
            self.writeSingleString(file, line)

            # Look up enthalpy of reaction
            # Define placeholder line (energy = 0); replace if necessary
            line = '\\\n' + len(startFirstLine) * ' ' + 'state.reactions[{0}].energy = 0;'.format(i)
            for j in range(len(enthalpiesList)):
                if reactionsList[i][0] == enthalpiesList[j][0]:
                    energy = Decimal(enthalpiesList[j][1]) * 1000 / Decimal(self.Avogadro)   # Convert from kJ/mol to J
                    energy = Decimal(energy) / Decimal(generalDict['heatcapacity'])          # Convert from J to K.dm3
                    line = '\\\n' + len(startFirstLine) * ' ' + 'state.reactions[{0}].energy = {1:.8e};'.format(i, energy)
                    break
            if i == (len(reactionsList) - 1):
                line = line + '}\n'
            else:
                line = line + '\\\n'
            self.writeSingleString(file, line)
 
        return True


    def RATES_UPDATE_BODY(self, file, generalDict, moleculesList, reactionsList, ratesList):
    # We're going to define how rate constants need to be updated

        particlecount = generalDict['particlecount']
        simulateheating = generalDict['simulateheating']
        concentration = self.calcConcentration(moleculesList)
        startFirstLine = '#define RATES_UPDATE_BODY {'
        body = startFirstLine
        updateCases = [] # Types of rates which require updating because of 
                         #   1) temperature dependence ('a', 'fd', 'ad', 'fp', 'ap', 'fpr', 'apr'), though only when SIMULATEHEATING is defined
                         #   2) conversion dependence ('fd', 'ad', 'fp', 'ap', 'fpr', 'apr')
                         #   3) dependence upon a second reaction process ('fpr', 'apr')
        
        # Populate updateCases
        if simulateheating != 0:
            updateCases.extend(['a', 'fd', 'ad', 'fp', 'ap', 'fpr', 'apr'])
        updateCases.extend(['fd', 'ad', 'fp', 'ap', 'fpr', 'apr'])
        updateCases.extend(['fpr', 'apr'])
        updateCases = list(set(updateCases)) # Remove duplicates
        
        # Go through the list of reactions
        for i in range(len(reactionsList)):

            # Look up the corresponding rate coefficient details
            for j in range(len(ratesList)):
                if reactionsList[i][1] == ratesList[j][0]:
                    rate = ratesList[j]
                    break

            # Describe a rate with Arrhenius parameters
            if rate[1] in updateCases:
                
                # Compile the string that describes the rate coefficient, append
                rcline = self.compileRateString(rate,particlecount,concentration,reactionsList[i],i)
                if body == startFirstLine:
                    body = body + rcline
                else:
                    body = body + '\\\n' + len(startFirstLine) * ' ' + rcline

            # When done iterating, terminate body with curly bracket, then write to output
            if i == (len(reactionsList) - 1):
                body = body + '}\n'
                self.writeSingleString(file, body)
 
        return True


    def MOLECULENAMES(self, file, moleculesList):
    # We're going to define each molecule name

        startFirstLine = '#define MOLECULENAMES {'

        # Go through the list of molecules, print each name
        for i in range(len(moleculesList)):
            
            # Look up the name
            name = moleculesList[i][1]
        
            # Print the name
            line = 'case {0}:\\\n'.format(name)
            if i == 0:
                line = startFirstLine + line
            else:
                line = len(startFirstLine) * ' ' + line
            self.writeSingleString(file, line)
            line = len(startFirstLine) * ' ' + '  ' + 'return ("{0}");'.format(name)
            if i == (len(moleculesList) - 1):
                line = line + '}\n'
            else:
                line = line + '\\\n'
            self.writeSingleString(file, line)

        return True


    def REACTIONNAMES(self, file, reactionsList):
    # We're going to define each reaction name

        startFirstLine = '#define REACTIONNAMES {'

        # Go through the list of molecules, print each name
        for i in range(len(reactionsList)):
            
            # Look up the name
            name = reactionsList[i][0]
        
            # Print the name
            line = 'case {0}:\\\n'.format(name)
            if i == 0:
                line = startFirstLine + line
            else:
                line = len(startFirstLine) * ' ' + line
            self.writeSingleString(file, line)
            line = len(startFirstLine) * ' ' + '  ' + 'return ("{0}");'.format(name)
            if i == (len(reactionsList) - 1):
                line = line + '}\n'
            else:
                line = line + '\\\n'
            self.writeSingleString(file, line)

        return True


    def DO_REACT_BODY(self, file, generalDict, moleculesList, reactionsList):
    # We're going to define the result of each reaction

        startFirstLine = '#define DO_REACT_BODY {'

        # Go through the list of reactions. For, each specify
        # - the case name
        # - the number of results (products)
        # - the name of the results (products)
        # - the arm lengths of the products ('poly' species have 1 arm, 'complex' species 2 ore more)
        for i in range(len(reactionsList)):

            # Determine the reaction type
            generalReaction = read.generalizeReaction(reactionsList[i], moleculesList)
            reactionHash = read.hashReaction(generalReaction)
            reactionType = read.reactionsDict[reactionHash]

            # Write the case name
            if i == 0:
                line = startFirstLine
            else:
                line = len(startFirstLine) * ' '
            line = line + 'case {0}:\\\n'.format(reactionsList[i][0])
            self.writeSingleString(file, line)

            # Write the number of results
            if reactionsList[i][5] != '':
                line = len(startFirstLine) * ' ' + '  no_of_res = 2;\\\n'
                self.writeSingleString(file, line)
            
            # Write the name of the first result
            p1 = read.isReactingComplex('',reactionsList[i][4],moleculesList)[1]
            line = len(startFirstLine) * ' ' + '  prod1_ind = {0};\\\n'.format(p1)
            self.writeSingleString(file, line)

            # Write the name of the second result
            if reactionsList[i][5] != '':
                p2 = read.isReactingComplex('',reactionsList[i][5],moleculesList)[1]
                line = len(startFirstLine) * ' ' + '  prod2_ind = {0};\\\n'.format(p2)
                self.writeSingleString(file, line)

            # Detail the arm length(s), etc. of the product(s). This part probably deserves its own function
            self.writeDetailedDO_REACT_BODY(file, generalDict, moleculesList, reactionsList[i], reactionType, startFirstLine)

            # Write a 'break'
            line = len(startFirstLine) * ' ' + '  break;'.format(p1)
            self.writeSingleString(file, line)
            if i == (len(reactionsList) - 1):
                line = '}\n'
            else:
                line = '\\\n'
            self.writeSingleString(file, line)

        return True


    def MWD_INITS(self, file, generalDict, moleculesList):
    # We're going to define the number of molecules for each reactant with a non-zero concentration

        startFirstLine = '#define MWD_INITS {'
        particlecount = generalDict['particlecount']
        concentration = self.calcConcentration(moleculesList)
        moleculesListNonzero = []

        for i in range(len(moleculesList)):
            if moleculesList[i][0] == 'simple' and moleculesList[i][2] != 0:
                moleculesListNonzero.append(moleculesList[i])

        for i in range(len(moleculesListNonzero)):
            count = int(round(moleculesListNonzero[i][2] * particlecount / concentration))
            if i == 0:
                line = startFirstLine + 'state.ms_cnts[{0}] = {1};'.format(moleculesListNonzero[i][1],count)
            else:
                line = len(startFirstLine) * ' ' + 'state.ms_cnts[{0}] = {1};'.format(moleculesListNonzero[i][1],count)
            self.writeSingleString(file, line)
            if i == (len(moleculesListNonzero) - 1):
                line = '}\n'
            else:
                line = '\\\n'
            self.writeSingleString(file, line)

        return True


    def writeSingleString(self, file, output):
    # A simple function that writes a single string to the output file

        try:
            f = open(file,'a')
            f.write(output)
            f.flush()
            f.close()
            return True
        except:
            print('\n An unknown error encountered while writing output. Aborting...')
            sys.exit()


    def calcConcentration(self, moleculesList):
    # Calculates the total concentration of reactive particles
        concentration = 0
        for i in range(len(moleculesList)):
            if moleculesList[i][0] == 'simple':
                concentration = float(Decimal(concentration) + Decimal(moleculesList[i][2]))
        return concentration


    def compileRateString(self, rate, particlecount, concentration, reaction, index):
        
        bimolCorStr = ''

        # Describe a rate with a fixed rate coefficient
        if rate[1] == 'f':
            
            if reaction[3] == '': # Unimolecular reaction, no rate adjustment needed
                rc = rate[2]
            else: # Bimolecular reaction, so we need to adjust the rate coefficient
                bimolCorStr = '((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * '
                rc = Decimal(rate[2]) * Decimal(concentration) / Decimal(particlecount)
                if reaction[2] == reaction[3]: # Bimolecular reaction with identical reactant, so additional multiplication by 2 needed
                    rc = rc * Decimal(2)
            rcline = 'state.reactions[{0}].rc = {1}{2:.8e};'.format(index, bimolCorStr, float(rc))

        # Describe a rate with Arrhenius parameters
        elif rate[1] == 'a':
            preexp = rate[2]
            exp = Decimal(-1) * Decimal(rate[3]) * Decimal(1000) / Decimal(self.gasconstant)
            if reaction[3] == '': # Unimolecular reaction, no rate adjustment needed
                pass
            else: # Bimolecular reaction, so we need to adjust the rate coefficient
                bimolCorStr = '((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * '
                preexp = Decimal(preexp) * Decimal(concentration) / Decimal(particlecount)
                if reaction[2] == reaction[3]: # Bimolecular reaction with identical reactant, so additional multiplication by 2 needed
                    preexp = preexp * Decimal(2)
            rcline = 'state.reactions[{0}].rc = {1}{2:.8e} * exp({3:.8e} / state.temp);'.format(index, bimolCorStr, float(preexp), float(exp))

        # Describe a rate with a fixed rate coefficient and damping
        elif rate[1] == 'fd':
            if reaction[3] == '': # Unimolecular reaction, no rate adjustment needed
                rc = rate[2]
            else: # Bimolecular reaction, so we need to adjust the rate coefficient
                bimolCorStr = '((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * '
                rc = Decimal(rate[2]) * Decimal(concentration) / Decimal(particlecount)
                if reaction[2] == reaction[3]: # Bimolecular reaction with identical reactant, so additional multiplication by 2 needed
                    rc = rc * Decimal(2)
            c1t0 = rate[3]
            c1t1 = rate[4]
            c2t0 = rate[5]
            c2t1 = rate[6]
            c3t0 = rate[7]
            c3t1 = rate[8]
            rcline = 'state.reactions[{0}].rc = {1}{2:.8e}'.format(index, bimolCorStr, float(rc))
            rcline = rcline + ' * exp(-1.0 * ( ({0:.4e}+{1:.4e}*state.temp)*state.conversion + ({2:.4e}+{3:.4e}*state.temp)*state.conversion*state.conversion + ({4:.4e}+{5:.4e}*state.temp)*state.conversion*state.conversion*state.conversion ) );'.format(c1t0,c1t1,c2t0,c2t1,c3t0,c3t1)

        # Describe a rate with Arrhenius parameters and damping
        elif rate[1] == 'ad':
            preexp = rate[2]
            exp = Decimal(-1) * Decimal(rate[3]) * Decimal(1000) / Decimal(self.gasconstant)
            if reaction[3] == '': # Unimolecular reaction, no rate adjustment needed
                pass
            else: # Bimolecular reaction, so we need to adjust the rate coefficient
                bimolCorStr = '((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * '
                preexp = Decimal(preexp) * Decimal(concentration) / Decimal(particlecount)
                if reaction[2] == reaction[3]: # Bimolecular reaction with identical reactant, so additional multiplication by 2 needed
                    preexp = preexp * Decimal(2)
            c1t0 = rate[4]
            c1t1 = rate[5]
            c2t0 = rate[6]
            c2t1 = rate[7]
            c3t0 = rate[8]
            c3t1 = rate[9]
            rcline = 'state.reactions[{0}].rc = {1}{2:.8e} * exp({3:.8e} / state.temp)'.format(index, bimolCorStr, float(preexp), float(exp))
            rcline = rcline + ' * exp(-1.0 * ( ({0:.4e}+{1:.4e}*state.temp)*state.conversion + ({2:.4e}+{3:.4e}*state.temp)*state.conversion*state.conversion + ({4:.4e}+{5:.4e}*state.temp)*state.conversion*state.conversion*state.conversion ) );'.format(c1t0,c1t1,c2t0,c2t1,c3t0,c3t1)

        # Describe a rate with a fixed rate coefficient and diffusion considered by a simple 'parallel' encounter pair model
        elif rate[1] == 'fp':
            # All parameters
            kc = rate[2]
            kd0 = rate[3]
            A = rate[4]
            B = rate[5]

            if reaction[3] == '': # Unimolecular reaction, no rate adjustment needed
                pass
            else: # Bimolecular reaction, so we need to adjust the rate coefficient
                bimolCorStr = '((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * '
                kc = Decimal(kc) * Decimal(concentration) / Decimal(particlecount)
                kd0 = Decimal(kd0) * Decimal(concentration) / Decimal(particlecount)
                if reaction[2] == reaction[3]: # Bimolecular reaction with identical reactant, so additional multiplication by 2 needed
                    kc = kc * Decimal(2)
                    kd0 = kd0 * Decimal(2)

            # Write lambda powers
            if A == 0:
                l0pow_line = '1'
                l1pow_line = '1'
            elif A == int(A):      # Speed up power calculation when possible
                l0pow_line = '((double)state.momentDist[0]'
                l1pow_line = '((double)state.momentDist[1]'
                for a in range(1, int(A)):
                    l0pow_line = l0pow_line + '*(double)state.momentDist[0]'
                    l1pow_line = l1pow_line + '*(double)state.momentDist[1]'
                l0pow_line = l0pow_line + ')'
                l1pow_line = l1pow_line + ')'
            else:
                l0pow_line = 'pow((double)state.momentDist[0], {0})'.format(A)
                l1pow_line = 'pow((double)state.momentDist[1], {0})'.format(A)

            num = '{0:.8e}*{1:.8e}*{2}'.format(float(kc), float(kd0), l0pow_line)
            denom_pt1 = '{0:.8e}*{1}*exp({2}/state.freeVolumeFraction)'.format(float(kc), l1pow_line, B)
            denom_pt2 = '{0:.8e}*{1}'.format(float(kd0), l0pow_line)
            rcline = 'state.reactions[{0}].rc = {1}({2}) / ( ({3}) + ({4}));'.format(index, bimolCorStr, num, denom_pt1, denom_pt2)

        # Describe a rate with a Arrhenius parameters and diffusion considered by a simple 'parallel' encounter pair model
        elif rate[1] == 'ap':
            # All parameters
            preexp = rate[2]
            exp = Decimal(rate[3]) * Decimal(1000) / Decimal(self.gasconstant)   # Caution: negative sign has been omitted as rewriting (1 / k = 1 / kc + 1 / kd) removes it
            kd0 = rate[4]
            A = rate[5]
            B = rate[6]

            if reaction[3] == '': # Unimolecular reaction, no rate adjustment needed
                pass
            else: # Bimolecular reaction, so we need to adjust the rate coefficient
                bimolCorStr = '((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * '
                preexp = Decimal(preexp) * Decimal(concentration) / Decimal(particlecount)
                kd0 = Decimal(kd0) * Decimal(concentration) / Decimal(particlecount)
                if reaction[2] == reaction[3]: # Bimolecular reaction with identical reactant, so additional multiplication by 2 needed
                    preexp = preexp * Decimal(2)
                    kd0 = kd0 * Decimal(2)

            # Write lambda powers
            if A == 0:
                l0pow_line = '1'
                l1pow_line = '1'
            elif A == int(A):      # Speed up power calculation when possible
                l0pow_line = '((double)state.momentDist[0]'
                l1pow_line = '((double)state.momentDist[1]'
                for a in range(1, int(A)):
                    l0pow_line = l0pow_line + '*(double)state.momentDist[0]'
                    l1pow_line = l1pow_line + '*(double)state.momentDist[1]'
                l0pow_line = l0pow_line + ')'
                l1pow_line = l1pow_line + ')'
            else:
                l0pow_line = 'pow((double)state.momentDist[0], {0})'.format(A)
                l1pow_line = 'pow((double)state.momentDist[1], {0})'.format(A)

            num = '{0:.8e}*{1:.8e}*{2}'.format(float(preexp), float(kd0), l0pow_line)
            denom_pt1 = '{0:.8e}*{1}*exp({2}/state.freeVolumeFraction)'.format(float(preexp), l1pow_line, B)
            denom_pt2 = '{0:.8e}*{1}*exp({2:.8e} / state.temp)'.format(float(kd0), l0pow_line, float(exp))
            rcline = 'state.reactions[{0}].rc = {1}({2}) / ( ({3}) + ({4}));'.format(index, bimolCorStr, num, denom_pt1, denom_pt2)

        # Describe a rate with a fixed rate coefficient and diffusion considered by a 'parallel' encounter pair model + 1 additional term
        elif rate[1] == 'fpr':
            # All parameters
            kc = rate[2]
            kd0 = rate[3]
            A = rate[4]
            B = rate[5]
            extRate = rate[6]
            species = rate[7]

            if reaction[3] == '': # Unimolecular reaction, no rate adjustment needed
                pass
            else: # Bimolecular reaction, so we need to adjust the rate coefficient
                bimolCorStr = '((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * '
                kc = Decimal(kc) * Decimal(concentration) / Decimal(particlecount)
                kd0 = Decimal(kd0) * Decimal(concentration) / Decimal(particlecount)
                extRate = Decimal(extRate) * Decimal(concentration) / Decimal(particlecount)
                if reaction[2] == reaction[3]: # Bimolecular reaction with identical reactant, so additional multiplication by 2 needed
                    kc = kc * Decimal(2)
                    extRate = kd0 * Decimal(2)
                    extRate = extRate * Decimal(2)

            # Write lambda powers
            if A == 0:
                rn = '1'
            elif A == int(A):      # Speed up power calculation when possible
                rn = '((double)state.momentDist[0]/(double)state.momentDist[1])'
                for a in range(1, int(A)):
                    rn = rn + '*((double)state.momentDist[0]/(double)state.momentDist[1])'
            else:
                rn = 'pow(((double)state.momentDist[0]/(double)state.momentDist[1]), {0})'.format(A)

            kd = '({0:.8e}*{1}*exp(-{2}/state.freeVolumeFraction))'.format(float(kd0), rn, B)
            kdr = '({0:.8e}*(state.ms_cnts[{1}]/AVOGADRO/state.volume))'.format(float(extRate), species)
            rcline = 'state.reactions[{0}].rc = {1}(1 / ((1/({2:.8e})) + (1/({3}+{4}))));'.format(index, bimolCorStr, kc, kd, kdr)

        # Describe a rate with a Arrhenius parameters and diffusion considered by a 'parallel' encounter pair model + 1 additional term
        elif rate[1] == 'apr':
            # All parameters
            preexp = rate[2]
            exp = Decimal(-1) * Decimal(rate[3]) * Decimal(1000) / Decimal(self.gasconstant)
            kd0 = rate[4]
            A = rate[5]
            B = rate[6]
            extRatePreexp = rate[7]
            extRateExp = Decimal(-1) * Decimal(rate[8]) * Decimal(1000) / Decimal(self.gasconstant)
            species = rate[9]

            if reaction[3] == '': # Unimolecular reaction, no rate adjustment needed
                pass
            else: # Bimolecular reaction, so we need to adjust the rate coefficient
                bimolCorStr = '((double)GLOBAL_MONOMER_PARTICLES/(double)state.localMonomerParticles) * '
                preexp = Decimal(preexp) * Decimal(concentration) / Decimal(particlecount)
                kd0 = Decimal(kd0) * Decimal(concentration) / Decimal(particlecount)
                extRatePreexp = Decimal(extRatePreexp) * Decimal(concentration) / Decimal(particlecount)
                if reaction[2] == reaction[3]: # Bimolecular reaction with identical reactant, so additional multiplication by 2 needed
                    preexp = preexp * Decimal(2)
                    kd0 = kd0 * Decimal(2)
                    extRatePreexp = extRatePreexp * Decimal(2)

            # Write lambda powers
            if A == 0:
                rn = '1'
            elif A == int(A):      # Speed up power calculation when possible
                rn = '((double)state.momentDist[0]/(double)state.momentDist[1])'
                for a in range(1, int(A)):
                    rn = rn + '*((double)state.momentDist[0]/(double)state.momentDist[1])'
            else:
                rn = 'pow(((double)state.momentDist[0]/(double)state.momentDist[1]), {0})'.format(A)

            kc = '({0:.8e} * exp({1:.8e} / state.temp))'.format(float(preexp),float(exp))
            kd = '({0:.8e}*{1}*exp(-{2}/state.freeVolumeFraction))'.format(float(kd0), rn, B)
            kdr = '({0:.8e} * exp({1:.8e} / state.temp)*(state.ms_cnts[{2}]/AVOGADRO/state.volume))'.format(float(extRatePreexp), float(extRateExp), species)
            rcline = 'state.reactions[{0}].rc = {1}(1 / ((1/({2})) + (1/({3}+{4}))));'.format(index, bimolCorStr, kc, kd, kdr)

        return rcline        


    def writeDetailedDO_REACT_BODY(self, file, generalDict, moleculesList, reaction, reactionType, startFirstLine):
    # Writes detailed DO_REACT_BODY information regarding how polymer lengths and complex arms lengths carry over
    # It is probably the least complicated to just write a case for each reaction type

        if reactionType == 'unimol_conversion':             # Simple   +  None     ->  Simple   +  None
            return True     # There are no lengths to define

        elif reactionType == 'decomposition':               # Simple   +  None     ->  Simple   +  Simple
            return True     # There are no lengths to define

        elif reactionType == 'combination_nonpoly':         # Simple   +  Simple   ->  Simple   +  None
            return True     # There are no lengths to define

        elif reactionType == 'transfer_nonpoly':            # Simple   +  Simple   ->  Simple   +  Simple
            return True     # There are no lengths to define

        elif reactionType == 'initiation2':                 # Simple   +  None     ->  Poly     +  None
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = 1;\\\n'
            self.writeSingleString(file, line)
            return True

        elif reactionType == 'initiation':                  # Simple   +  Simple   ->  Poly     +  None
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = 1;\\\n'
            self.writeSingleString(file, line)
            return True

        elif reactionType == 'polymeric_change':            #  Poly     +  None     ->  Poly     +  None
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = react1_lens[0];\\\n'
            self.writeSingleString(file, line)
            return True

        elif reactionType == 'depropagation':               # Poly     +  None     ->  Poly     +  Simple
            if reaction[5] in generalDict['monomernames']:
                line = len(startFirstLine) * ' ' + '  prod1_lens[0] = react1_lens[0] - 1;\\\n'    # Depropagation
            else:
                line = len(startFirstLine) * ' ' + '  prod1_lens[0] = react1_lens[0];\\\n'        # Cleavage of a non-monomeric species
            self.writeSingleString(file, line)
            return True

        elif reactionType == 'propagation':                 # Poly     +  Simple   ->  Poly     +  None
            if reaction[3] in generalDict['monomernames']:
                line = len(startFirstLine) * ' ' + '  prod1_lens[0] = react1_lens[0] + 1;\\\n'    # Propagation
            else:
                line = len(startFirstLine) * ' ' + '  prod1_lens[0] = react1_lens[0];\\\n'        # Addition of a non-monomeric species
            self.writeSingleString(file, line)
            return True

        elif reactionType == 'chaintransfer':               # Poly     +  Simple   ->  Poly     +  Simple
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = react1_lens[0];\\\n'
            self.writeSingleString(file, line)
            return True

        elif reactionType == 'transfer_with_poly_creation':   # Poly     +  Simple   ->  Poly     +  Poly
                                                              # This is a tough case. To the user, it is not clear that second product is seen as the newly 
                                                              # initiated polymer. This is to be either explicitely communicated (in documentation
                                                              # and a warning in this program) or improved upon in a future release (through a revised 
                                                              # input scheme).
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = react1_lens[0];\\\n'
            self.writeSingleString(file, line)
            line = len(startFirstLine) * ' ' + '  prod2_lens[0] = 1;\\\n'
            self.writeSingleString(file, line)
            return True

        elif reactionType == 'combination':                 # Poly     +  Poly     ->  Poly     +  None
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = react1_lens[0] + react2_lens[0];\\\n'
            self.writeSingleString(file, line)
            return True

        elif reactionType == 'disproportionation':          # Poly     +  Poly     ->  Poly     +  Poly
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = react1_lens[0];\\\n'
            self.writeSingleString(file, line)
            line = len(startFirstLine) * ' ' + '  prod2_lens[0] = react2_lens[0];\\\n'
            self.writeSingleString(file, line)
            return True

        elif reactionType == '2polys_to_1complex':          # Poly     +  Poly     ->  Complex  +  None
            # No need to look up the arm number -- a newly formed complex always has 2 arms
            line = len(startFirstLine) * ' ' + '  prod1_arms = 2;\\\n'
            self.writeSingleString(file, line)
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = react1_lens[0];\\\n'
            self.writeSingleString(file, line)
            line = len(startFirstLine) * ' ' + '  prod1_lens[1] = react2_lens[0];\\\n'
            self.writeSingleString(file, line)
            return True
            
        elif reactionType == 'consolidate_type1':           # Complex  +  None     ->  Poly     +  None
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = '
            self.writeSingleString(file, line)

            # Look up the arm number of reacting complex species
            r1 = read.isReactingComplex('',reaction[2],moleculesList)[1]
            for i in range(len(moleculesList)):
                if r1 == moleculesList[i][1]:
                    r1_arms = moleculesList[i][2]

            # The product length is the sum of the arm lengths            
            for i in range(r1_arms):
                line = 'react1_lens[{0}]'.format(i)
                if i != (r1_arms - 1):
                    line = line + ' + '
                else:
                    line = line + ';\\\n'
                self.writeSingleString(file, line)
            return True

        elif reactionType == 'complex_cleave_2polys':       # Complex  +  None     ->  Poly     +  Poly
                                                            # This is another tough case: which arm becomes which product? As currently coded, the 
                                                            # the index becomes the first polymer species. This is not clear to the user though, so it 
                                                            # is probably best that this is either explicitely communicated (in documentation
                                                            # and a warning in this program) or improved upon in a future release (through a revised 
                                                            # input scheme).
            # No need to look up the arm number -- a complex broken into two polys always has 2 arms
            # We do need to look up the arm index though
            r1_index = read.isReactingComplex('',reaction[2],moleculesList)[2]   # Can be '0' or '1'
            
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = react1_lens[{0}];\\\n'.format(r1_index)
            self.writeSingleString(file, line)
            line = len(startFirstLine) * ' ' + '  prod2_lens[0] = react1_lens[{0}];\\\n'.format(abs(r1_index - 1))
            self.writeSingleString(file, line)
            return True

        elif reactionType == 'complex_change':              # Complex  +  None     ->  Complex  +  None
                                                            # Does this type of reaction even matter (will it be used)?
            # Look up the arm number of product (which equals that of the reactant)
            p1 = read.isReactingComplex('',reaction[4],moleculesList)[1]
            for i in range(len(moleculesList)):
                if p1 == moleculesList[i][1]:
                    p1_arms = moleculesList[i][2]

            line = len(startFirstLine) * ' ' + '  prod1_arms = {0};\\\n'.format(p1_arms)
            self.writeSingleString(file, line)
            for i in range(p1_arms):
                line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react1_lens[{0}];\\\n'.format(i)
                self.writeSingleString(file, line)
            return True

        elif reactionType == 'depropagation_complex':       # Complex  +  None     ->  Complex  +  Simple
            # Look up the arm number of product (which equals that of the reactant)
            p1 = read.isReactingComplex('',reaction[4],moleculesList)[1]
            for i in range(len(moleculesList)):
                if p1 == moleculesList[i][1]:
                    p1_arms = moleculesList[i][2]

            # Look up the index of the depropagating arm
            r1_index = read.isReactingComplex('',reaction[2],moleculesList)[2]
            
            line = len(startFirstLine) * ' ' + '  prod1_arms = {0};\\\n'.format(p1_arms)
            self.writeSingleString(file, line)
            for i in range(p1_arms):
                if i == r1_index:
                    if reaction[5] in generalDict['monomernames']:
                        line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react1_lens[{0}] - 1;\\\n'.format(i)    # Depropagation
                    else:
                        line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react1_lens[{0}];\\\n'.format(i)        # Cleavage of a non-monomeric species
                else:
                    line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react1_lens[{0}];\\\n'.format(i)
                self.writeSingleString(file, line)
            return True
            
        elif reactionType == 'complex_cleave_1poly':        # Complex  +  None     ->  Complex  +  Poly
            # Look up the arm number of product (which equals that of the reactant minus 1)
            p1 = read.isReactingComplex('',reaction[4],moleculesList)[1]
            for i in range(len(moleculesList)):
                if p1 == moleculesList[i][1]:
                    p1_arms = moleculesList[i][2]
            line = len(startFirstLine) * ' ' + '  prod1_arms = {0};\\\n'.format(p1_arms)
            self.writeSingleString(file, line)

            # Look up the index of the arm to cleave
            r1_index = read.isReactingComplex('',reaction[2],moleculesList)[2]

            # Define each reactant complex arm as a product complex arm, save for the arm to cleave
            indexCounter = 0
            for i in range(p1_arms):
                if i == r1_index:
                    indexCounter = indexCounter + 1    # Skip the index of the arm to cleave
                line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react1_lens[{1}];\\\n'.format(i, indexCounter)
                self.writeSingleString(file, line)        
                indexCounter = indexCounter + 1

            # The cleaved arm becomes a poly
            line = len(startFirstLine) * ' ' + '  prod2_lens[0] = react1_lens[{0}];\\\n'.format(r1_index)
            self.writeSingleString(file, line)
            return True

        elif reactionType == 'propagation_complex':         # Complex  +  Simple   ->  Complex  +  None
            # Look up the arm number of product (which equals that of the reactant)
            p1 = read.isReactingComplex('',reaction[4],moleculesList)[1]
            for i in range(len(moleculesList)):
                if p1 == moleculesList[i][1]:
                    p1_arms = moleculesList[i][2]

            # Look up the index of the propagating arm
            r1_index = read.isReactingComplex('',reaction[2],moleculesList)[2]
            
            line = len(startFirstLine) * ' ' + '  prod1_arms = {0};\\\n'.format(p1_arms)
            self.writeSingleString(file, line)
            for i in range(p1_arms):
                if i == r1_index:
                    if reaction[3] in generalDict['monomernames']:
                        line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react1_lens[{0}] + 1;\\\n'.format(i)    # Propagation
                    else:
                        line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react1_lens[{0}];\\\n'.format(i)        # Cleavage of a non-monomeric species
                else:
                    line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react1_lens[{0}];\\\n'.format(i)
                self.writeSingleString(file, line)
            return True

        elif reactionType == 'consolidate_type3':           # Complex  +  Poly     ->  Poly     +  None
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = '
            self.writeSingleString(file, line)

            # Look up the arm number of reacting complex species
            r1 = read.isReactingComplex('',reaction[2],moleculesList)[1]
            for i in range(len(moleculesList)):
                if r1 == moleculesList[i][1]:
                    r1_arms = moleculesList[i][2]

            # The product length is the sum of the arm lengths, plus the arm length of the poly species
            for i in range(r1_arms):
                line = 'react1_lens[{0}]'.format(i)
                if i != (r1_arms - 1):
                    line = line + ' + '
                else:
                    line = line + ' + react2_lens[0];\\\n'
                self.writeSingleString(file, line)
            return True

        elif reactionType == 'complex_expand_1':            # Complex  +  Poly     ->  Complex  +  None
            # Look up the arm number of product (which equals that of the reactant plus 1)
            p1 = read.isReactingComplex('',reaction[4],moleculesList)[1]
            for i in range(len(moleculesList)):
                if p1 == moleculesList[i][1]:
                    p1_arms = moleculesList[i][2]
            line = len(startFirstLine) * ' ' + '  prod1_arms = {0};\\\n'.format(p1_arms)
            self.writeSingleString(file, line)

            # Look up the index where the new arm should be placed on
            r1_index = read.isReactingComplex('',reaction[2],moleculesList)[2]

            # Define each reactant complex arm as a product complex arm while inserting the new arm on the correct index
            indexCounter = 0
            for i in range(p1_arms):
                if i == r1_index:   # Insert the new arm
                    line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react2_lens[0];\\\n'.format(i)
                    self.writeSingleString(file, line)
                else:                    
                    line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react1_lens[{1}];\\\n'.format(i, indexCounter)
                    self.writeSingleString(file, line)
                    indexCounter = indexCounter + 1

            return True

        elif reactionType == 'consolidate_type2':           # Complex  +  Complex  ->  Poly     +  None
            line = len(startFirstLine) * ' ' + '  prod1_lens[0] = '
            self.writeSingleString(file, line)

            # Look up the arm numbers for two reacting complex species
            r1 = read.isReactingComplex('',reaction[2],moleculesList)[1]
            for i in range(len(moleculesList)):
                if r1 == moleculesList[i][1]:
                    r1_arms = moleculesList[i][2]
            r2 = read.isReactingComplex('',reaction[3],moleculesList)[1]
            for i in range(len(moleculesList)):
                if r2 == moleculesList[i][1]:
                    r2_arms = moleculesList[i][2]

            # The product length is the sum of the two arm lengths
            for i in range(r1_arms):
                line = 'react1_lens[{0}] + '.format(i)
                self.writeSingleString(file, line)
            for i in range(r2_arms):
                line = 'react2_lens[{0}]'.format(i)
                if i != (r2_arms - 1):
                    line = line + ' + '
                else:
                    line = line + ';\\\n'
                self.writeSingleString(file, line)
            return True

        elif reactionType == 'complex_combination':         # Complex  +  Complex  ->  Complex  +  None
            # Look up the arm numbers for two reacting complex species
            r1 = read.isReactingComplex('',reaction[2],moleculesList)[1]
            for i in range(len(moleculesList)):
                if r1 == moleculesList[i][1]:
                    r1_arms = moleculesList[i][2]
            r2 = read.isReactingComplex('',reaction[3],moleculesList)[1]
            for i in range(len(moleculesList)):
                if r2 == moleculesList[i][1]:
                    r2_arms = moleculesList[i][2]

            # The arms of the two reactant complexes are combined in the product complex
            for i in range(r1):
                line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react1_lens[{0}];\\\n'.format(i)
                self.writeSingleString(file, line)
            for i in range(r1,r1+r2):
                line = len(startFirstLine) * ' ' + '  prod1_lens[{0}] = react2_lens[{1}];\\\n'.format(i, i - r1)
                self.writeSingleString(file, line)
            return True
