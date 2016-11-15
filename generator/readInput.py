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

import linecache
import sys
import re

# For unit testing
try:
    from StringIO import StringIO  # For Python 2.x
except ImportError:
    from io import StringIO  # For Python 3.x

class readInput(object):
    
    def __init__(self):
        self.groups = ['general','molecules','reactions','ratecoefficients']
        # Reaction hash explained:
            # None = 1, Simple = 2, Poly = 3, Complex = 4
            # Reactant1: x**11, reactant2: x**7, product1: x**5, product2: x**3;
            # Hash is the sum of the four values 
        self.reactionsDict = {
        2082:    'unimol_conversion',              #  Simple   +  None     ->  Simple   +  None
        2089:    'decomposition',                  #  Simple   +  None     ->  Simple   +  Simple
        2209:    'combination_nonpoly',            #  Simple   +  Simple   ->  Simple   +  None
        2216:    'transfer_nonpoly',               #  Simple   +  Simple   ->  Simple   +  Simple
        2293:    'initiation2',                    #  Simple   +  None     ->  Poly     +  None
        2420:    'initiation',                     #  Simple   +  Simple   ->  Poly     +  None
        177392:  'polymeric_change',               #  Poly     +  None     ->  Poly     +  None
        177399:  'depropagation',                  #  Poly     +  None     ->  Poly     +  Simple   #  Whether this is depropagation (monomer cleavage) or other cleavage will be determined later
        177519:  'propagation',                    #  Poly     +  Simple   ->  Poly     +  None     #  Whether this is propagation (monomer addition) or other addition will be determined later
        177526:  'chaintransfer',                  #  Poly     +  Simple   ->  Poly     +  Simple
        177545:  'transfer_with_poly_creation',    #  Poly     +  Simple   ->  Poly     +  Poly
        179578:  'combination',                    #  Poly     +  Poly     ->  Poly     +  None
        179604:  'disproportionation',             #  Poly     +  Poly     ->  Poly     +  Poly
        180359:  '2polys_to_1complex',             #  Poly     +  Poly     ->  Complex  +  None
        4194549: 'consolidate_type1',              #  Complex  +  None     ->  Poly     +  None
        4194575: 'complex_cleave_2polys',          #  Complex  +  None     ->  Poly     +  Poly
        4195330: 'complex_change',                 #  Complex  +  None     ->  Complex  +  None
        4195337: 'depropagation_complex',          #  Complex  +  None     ->  Complex  +  Simple   #  Whether this is depropagation (monomer cleavage) or other cleavage will be determined later
        4195356: 'complex_cleave_1poly',           #  Complex  +  None     ->  Complex  +  Poly
        4195457: 'propagation_complex',            #  Complex  +  Simple   ->  Complex  +  None     #  Whether this is depropagation (monomer cleavage) or other cleavage will be determined later
        4196735: 'consolidate_type3',              #  Complex  +  Poly     ->  Poly     +  None
        4197516: 'complex_expand_1',               #  Complex  +  Poly     ->  Complex  +  None
        4210932: 'consolidate_type2',              #  Complex  +  Complex  ->  Poly     +  None
        4211713: 'complex_combination',            #  Complex  +  Complex  ->  Complex  +  None
        }


    def readGeneral(self, file):
        # Perhaps this function should be chopped up into smaller ones        

        # Initialise defaults
        generalDict = {
        'particlecount':     -1,         # Int, required, no default exists
        'monomernames':      'defnames', # String, required, no default exists
        'seed':              -1,         # Int, not required, no default exists
        'calcdist':          0,          # Int, not required, default is 0 (= off)
        
        'syncsimtime':       0,          # Int, required if 'syncevents' is not specified, default is 0 (milliseconds)
        'syncevents':        0,          # Int, required if 'syncsimtime' is not specified, default is 0 (events)
        
        'maxsimtime':        0,          # Int, not required, default is 0 (seconds)
        'maxwalltime':       0,          # Int, not required, default is 0 (minutes)
        'maxevents':         0,          # Int, not required, default is 0 (events)
        'maxconversion':     0.0,        # Float, not required, default is 0.0 (conversion)
        
        'temperature':       0.0,        # Float, not required unless Arrhenius parameters are specified or free volume based diffusion control is enabled (checks performed elsewhere) (Kelvin)
        'starttemperature':  0.0,        # Float, not required, default is 0.0 (Kelvin)
        
        'simulateheating':   0,          # Int, not required, default is 0 (= off)
        'coolingrate':       0.0,        # Float, not required, default is 0.0 (s^-1)
        'heatcapacity':      0.0,        # Float, not required, default is 0.0 (J / (dm^3 K))
        
        'freevolume':        0,          # Int, not required, default is 0 (= off)
        'recalcconversion':  0,          # Int, not required, default is 0 (= off)
        'longchainsupport':  0,          # Int, not required, default is 0 (= off)
        'usefactor2':        1           # Int, not required, default is 1 (= on)
        }
        reqLength = [2]
        reqParams = ['particlecount','monomernames','maxwalltime','syncsimtime','syncevents']

        # Read lines of interest, store the parameters they contain in list
        inputParams = []
        linesList = self.readGroup(file,'general')
        if linesList == None:
            self.missingGroupError('general')
            return -1
        elif linesList == -1:
            return -1
        for i in linesList:
            line = linecache.getline(file, i).strip(' \t\n\r').split()
            line[0] = line[0].lower()
            # Check input validity
            if len(line) not in reqLength:
                self.lineError(linecache.getline(file, i),'Invalid number of arguments')
                return -1
            if line[0] not in generalDict:
                self.lineError(linecache.getline(file, i),('Unknown parameter "{0}"'.format(line[0])))
                return -1
            # Append to list of parameters that have bee input. In case of required input, strike off list of required inputs
            inputParams.append(line)
            if line[0] in reqParams:
                if line[0] == 'syncsimtime':
                    if 'syncevents' in reqParams:
                        reqParams.remove('syncevents')
                if line[0] == 'syncevents':
                    if 'syncsimtime' in reqParams:
                        reqParams.remove('syncsimtime')
                reqParams.remove(line[0])

        # Update dictionary with parameters read in
        for i in range(len(inputParams)):
            # Convert input value to right type (i.e., the type of the default value), abort when a wrong type is provided
            if isinstance(generalDict[inputParams[i][0]], float) == True:
                try:
                    inputParams[i][1] = float(inputParams[i][1])
                except ValueError:
                    self.valueTypeError('parameter',inputParams[i][0],inputParams[i][1],'Expected a float (although an int is allowed)')
                    return -1
            elif isinstance(generalDict[inputParams[i][0]], int) == True:
                try:
                    inputParams[i][1] = int(inputParams[i][1])
                except ValueError:  # If we cannot convert to int directly, try converting to float first
                    try:
                        inputParams[i][1] = int(float(inputParams[i][1]))
                    except ValueError:
                        self.valueTypeError('parameter',inputParams[i][0],inputParams[i][1],'Expected an int (although a float is allowed; it will be converted to int)')
                        return -1
            # Update dictionary
            generalDict[inputParams[i][0]] = inputParams[i][1]
        
        # Check if all required input parameters have been specified
        if len(reqParams) != 0:
            self.missingInputError('general',reqParams)
            return -1

        # Convert monomers to list, check if it is a list
        generalDict['monomernames'] = generalDict['monomernames'].split(',')
        if isinstance(generalDict['monomernames'],list) == False:
            print('\nError interpreting "monomername"'.format(seedMinVal,seedMaxVal))
            return -1

     # Check if the specified input is allowed
        # Check if the seed isn't too small or large
        seedMinVal = 0
        seedMaxVal = 2**32 - 1
        if generalDict['seed'] != -1 and (generalDict['seed'] < seedMinVal or generalDict['seed'] > seedMaxVal):
            print('\nThe specified seed is outside of the allowed range of {0} - {1}'.format(seedMinVal,seedMaxVal))
            return -1

        # Check if at least one stopping criterion has been provided
        if (generalDict['maxsimtime'] == 0) and (generalDict['maxwalltime'] == 0) and (generalDict['maxevents'] == 0) and (generalDict['maxconversion'] == 0):
            print('\nNo ending criterion has been provided. Please provide one or more')
            return -1

        # Check if maxsimtime and syncsimtime have a valid values
        if (generalDict['maxsimtime'] != 0):
            if (generalDict['maxsimtime'] < 0):
                print('\nInvalid value specified for the parameter "maxsimtime"'.format(seedMinVal,seedMaxVal))
                return -1
        if (generalDict['syncsimtime'] < 0):
            print('\nInvalid value specified for the parameter "syncsimtime"')
            return -1

        # Check if maxwalltime has a valid value
        if (generalDict['maxwalltime'] != 0):
            if (generalDict['maxwalltime'] < 0):
                print('\nInvalid value specified for the parameter "maxwalltime"')
                return -1

        # Check if maxevents and syncevents have a valid values
        if (generalDict['maxevents'] != 0):
            if (generalDict['maxevents'] < 0):
                print('\nInvalid value specified for the parameter "maxevents"')
                return -1
        if (generalDict['syncevents'] < 0):
            print('\nInvalid value specified for the parameter "syncsimtime"')
            return -1

        # Check if conversion has a valid value
        if (generalDict['maxconversion'] != 0):
            if (generalDict['maxconversion'] < 0) or (generalDict['maxconversion'] >= 1):
                print('\nInvalid value specified for the parameter "maxconversion"')
                return -1

        # Check whether the temperature is not negative
        if generalDict['temperature'] < 0.0:
            print('\nInvalid temperature "{0}" specified'.format(generalDict['temperature']))
            return -1

        # Check whether the starttemperature is not negative
        if generalDict['starttemperature'] < 0.0:
            print('\nInvalid starttemperature "{0}" specified'.format(generalDict['starttemperature']))
            return -1
        elif generalDict['starttemperature'] != generalDict['temperature'] and generalDict['coolingrate'] == 0.0:
            print('\nWarning: the start temperature differs from base temperature, but no cooling rate was specified')

        # Check if 'simulateheating' has a valid value
        if generalDict['simulateheating'] not in [0,1]:
            print('\nInvalid value specified for "simulateheating"')
            return -1

        # Check whether the heat capacity has a valid value
        if generalDict['heatcapacity'] < 0.0:
            print('\nInvalid heatcapacity "{0}" specified'.format(generalDict['heatcapacity']))
            return -1
        elif generalDict['heatcapacity'] == 0.0 and generalDict['simulateheating'] != 0:
            print('\nSimulation of heating has been requested, but no heat capacity has been specified')
            return -1

        # Check whether the cooling rate has a valid value
        if generalDict['coolingrate'] < 0.0:
            print('\nInvalid coolingrate "{0}" specified'.format(generalDict['coolingrate']))
            return -1

        # Check if 'calcdist' has a valid value
        if generalDict['calcdist'] not in [0,1]:
            print('\nInvalid value specified for "calcdist"')
            return -1

        # Check if 'freevolume' has a valid value
        if generalDict['freevolume'] not in [0,1]:
            print('\nInvalid value specified for "freevolume"')
            return -1
        
        # Check if 'recalcconversion' has a valid value
        if generalDict['recalcconversion'] not in [0,1]:
            print('\nInvalid value specified for "recalcconversion"')
            return -1
        
        # Check if 'longchainsupport' has a valid value
        if generalDict['longchainsupport'] not in [0,1]:
            print('\nInvalid value specified for "longchainsupport"')
            return -1
        
        # Check if 'longchainsupport' has a valid value
        if generalDict['usefactor2'] not in [0,1]:
            print('\nInvalid value specified for "usefactor2"')
            return -1

        return generalDict


    def readMolecules(self, file, generalDict):

        # Initialise defaults
        reqLengthSimple = [2,3]
        reqLengthPoly = [2]
        reqLengthComplex = [3]
        specifiedmonomers = generalDict['monomernames']
        monomer = []
        simple = []
        poly = []
        complex = []

        # Read lines of interest, store each molecule in the correct list
        linesList = self.readGroup(file,'molecules')
        if linesList == None:
            self.missingGroupError('molecules')
            return -1
        elif linesList == -1:
            return -1
        for i in linesList:
            line = linecache.getline(file, i).strip(' \t\n\r').split()
            line[0] = line[0].lower()

            # Check input validity
            if line[0] not in ['simple','poly','complex']:
                self.lineError(linecache.getline(file, i),'Unknown molecule type "{0}"'.format(line[0]))
                return -1
            if (line[0] == 'simple') and (len(line) not in reqLengthSimple):
                self.lineError(linecache.getline(file, i),'Invalid number of arguments')
                return -1
            if line[0] == 'poly' and (len(line) not in reqLengthPoly):
                self.lineError(linecache.getline(file, i),'Invalid number of arguments')
                return -1
            if line[0] == 'complex' and (len(line) not in reqLengthComplex):
                self.lineError(linecache.getline(file, i),'Invalid number of arguments')
                return -1

            # Initialise 'simple' concentrations, 'complex' arm lengths
            if line[0] == 'simple' and len(line) == reqLengthSimple[0]:
                line.append(0.0)
            elif line[0] == 'simple' and len(line) == reqLengthSimple[1]:
                try:
                    line[2] = float(line[2])
                except ValueError:
                    self.lineError(linecache.getline(file, i),'Invalid concentration specified')
                    return -1
            elif line[0] == 'complex':
                try:
                    line[2] = int(line[2])
                except ValueError:
                    self.lineError(linecache.getline(file, i),'Invalid number of arms specified')
                    return -1

            # Append to the correct list, but first check for duplicates
            combined = monomer + simple + poly + complex
            for j in range(len(combined)):
                if line[1] == combined[j][1]:
                    self.duplicateMolError(line[1])
                    return -1
            if line[0] == 'simple':
                if line[1] in specifiedmonomers:
                    monomer.append(line)
                else:
                    simple.append(line)
            elif line[0] == 'poly':
                poly.append(line)
            elif line[0] == 'complex':
                complex.append(line)

        # Check if all monomers named in 'general' have been found in 'molecules'
        for i in range(len(specifiedmonomers)):
            for j in range(len(monomer)):
                if specifiedmonomers[i] == monomer[j][1]:
                    break
                if j == (len(monomer) - 1):
                    print('\n"{0}" was specified as a monomer but has not been defined in the group "molecules"'.format(specifiedmonomers[i]))
                    return -1

        return (monomer + simple + poly + complex)


    def readReactions(self, file, moleculesList, ratesList):
        # Initialise defaults
        reqLength = [6, 8, 10]
        reactions = []
        warningsDisplayed = []  # Each reactionType-warning will only be have to be shown once

        # Read lines of interest, store each reaction in the list
        linesList = self.readGroup(file,'reactions')
        if linesList == None:
            self.missingGroupError('reactions')
            return -1
        elif linesList == -1:
            return -1
        for i in linesList:
            line = linecache.getline(file, i).strip(' \t\n\r').split()

            # Let's start with checking the length, and whether the third item is the '=' character
            if (len(line) not in reqLength) or (line[2] != '='):
                self.lineError(linecache.getline(file, i),'Invalid reaction structure')
                return -1

            # Check whether the rate coefficient is known
            for j in range(len(ratesList)):
                if line[1] == ratesList[j][0]:
                    break
                elif j == (len(ratesList) - 1):
                    self.lineError(linecache.getline(file, i),'Unknown rate coefficient')
                    return -1
            
            # Parse interconversions
            if len(line) == 6:

                # Additional checks regarding reaction syntax and molecule names
                if line[4] != '->':
                    self.lineError(linecache.getline(file, i),'Invalid reaction syntax')
                    return -1
                ans = self.checkSpecies(linecache.getline(file, i), line[3], 'r', moleculesList) # Check reactant1
                if ans == -1:
                    return  -1
                ans = self.checkSpecies(linecache.getline(file, i), line[5], 'p', moleculesList) # Check product1
                if ans == -1:
                    return  -1

                # Get rid of unneeded strings and introduce an empty string for the second reactant/product
                line.insert(4,'')       # Insert '' for reactant2
                line.append('')         # Insert '' for product2
                line.remove('=')        # Delete '='
                line.remove('->')       # Delete '->'

            # Parse combinations and decompositions
            elif (len(line) == 8):

                # Additional checks regarding reaction syntax and molecule names
                if (line[4] not in ['+','->']) or (line[6] not in ['+','->']) or (line[4] == line[6]):
                    self.lineError(linecache.getline(file, i),'Invalid reaction syntax')
                    return -1
                ans = self.checkSpecies(linecache.getline(file, i), line[3], 'r', moleculesList) # Check reactant1
                if ans == -1:
                    return  -1
                role = ('r' if line[4] == '+' else 'p')
                ans = self.checkSpecies(linecache.getline(file, i), line[5], role, moleculesList) # Check reactant2/product1
                if ans == -1:
                    return  -1
                ans = self.checkSpecies(linecache.getline(file, i), line[7], 'p', moleculesList) # Check product1/product2
                if ans == -1:
                    return  -1
                
                # Get rid of unneeded strings and introduce an empty string for the second reactant/product
                if line[4] == '+':
                    line.append('')     # Insert '' for product2
                if line[4] == '->':
                    line.insert(4,'')   # Insert '' for reactant2
                line.remove('+')        # Delete '='
                line.remove('=')        # Delete '='
                line.remove('->')       # Delete '->'

            # Parse transfers
            elif (len(line) == 10):

                # Additional checks regarding reaction syntax and molecule names
                if (line[4] != '+') or (line[6] != '->') or (line[8] != '+'):
                    self.lineError(linecache.getline(file, i),'Invalid reaction syntax')
                    return -1
                ans = self.checkSpecies(linecache.getline(file, i), line[3], 'r', moleculesList) # Check reactant1
                if ans == -1:
                    return  -1
                ans = self.checkSpecies(linecache.getline(file, i), line[5], 'r', moleculesList) # Check reactant2
                if ans == -1:
                    return  -1
                ans = self.checkSpecies(linecache.getline(file, i), line[7], 'p', moleculesList) # Check product1
                if ans == -1:
                    return  -1
                ans = self.checkSpecies(linecache.getline(file, i), line[9], 'p', moleculesList) # Check product2
                if ans == -1:
                    return  -1
                
                # Get rid of unneeded strings and introduce an empty string for the second reactant/product
                line.remove('+')    # Delete first '+'
                line.remove('+')    # Delete second '+'
                line.remove('=')    # Delete '='
                line.remove('->')   # Delete '->'

            # Generalize reaction, then order reactants and products according to complex -> poly -> simple -> none
            generalReaction = self.generalizeReaction(line, moleculesList)
            if generalReaction == -1:
                return -1
            if generalReaction[0] < generalReaction[1]:
                generalReaction[0], generalReaction[1] = generalReaction[1], generalReaction[0]
                line[2], line[3] = line[3], line[2]
            if generalReaction[2] < generalReaction[3]:
                generalReaction[2], generalReaction[3] = generalReaction[3], generalReaction[2]
                line[4], line[5] = line[5], line[4]

            # Hash the reaction and check if it exists in the dictionary
            reactionHash = self.hashReaction(generalReaction)
            if reactionHash not in self.reactionsDict:
                translation = ['None','Simple','Poly','Complex']
                problem = 'Reactions of the type "{0} + {1} -> {2} + {3}" are not supported'.format(translation[generalReaction[0]-1],translation[generalReaction[1]-1],translation[generalReaction[2]-1],translation[generalReaction[3]-1])
                self.reactionTypeError(linecache.getline(file, i), problem)
                return -1

            # For reactions involving complexes, check if the indices and arm numbers make sense
            # Request arm lengths and indices
            r1_name = self.isReactingComplex('',line[2],moleculesList)[1]
            if r1_name == -1:
                return -1
            if generalReaction[0] == 4:
                for j in range(len(moleculesList)):
                    if r1_name == moleculesList[j][1]:
                        r1_arms = moleculesList[j][2]
                        break
                r1_index = self.isReactingComplex('',line[2],moleculesList)[2]
                if r1_index == -1:
                    return -1
            r2_name = self.isReactingComplex('',line[3],moleculesList)[1]
            if r2_name == -1:
                return -1
            if generalReaction[1] == 4:
                for j in range(len(moleculesList)):
                    if r2_name == moleculesList[j][1]:
                        r2_arms = moleculesList[j][2]
                        break
                r2_index = self.isReactingComplex('',line[3],moleculesList)[2]
                if r2_index == -1:
                    return -1
            p1_name = line[4]
            if generalReaction[2] == 4:
                for j in range(len(moleculesList)):
                    if p1_name == moleculesList[j][1]:
                        p1_arms = moleculesList[j][2]
                        break
            p2_name = line[5]
            if generalReaction[3] == 4:
                for j in range(len(moleculesList)):
                    if p2_name == moleculesList[j][1]:
                        p2_arms = moleculesList[j][2]
                        break

            # Perform a check for each reaction type that can be a problem
            reactionType = self.reactionsDict[reactionHash]            
            if reactionType == 'transfer_with_poly_creation':  # Print a warning
                if reactionType not in warningsDisplayed:
                    print('\nWarning: a reaction of the type "Poly  +  Simple  ->  Poly  +  Poly" has been encountered. Although Simply is able to simulate such reactions without any problems, it cannot tell which of the two product polymer species is which. By design, it assumes the first product is the terminated chain and the second product is newly initiated chain. Please set up your input accordingly; swap the two product species if necessary.')
                    warningsDisplayed.append(reactionType)
            elif reactionType == '2polys_to_1complex':
                if p1_arms != 2:
                    self.reactionTypeError(linecache.getline(file, i), 'Product complex is specified as having more than 2 arms')
                    return -1
            elif reactionType == 'complex_cleave_2polys':
                if r1_arms != 2:
                    self.reactionTypeError(linecache.getline(file, i), 'Reactant complex is specified as having more than 2 arms')
                    return -1
                if reactionType not in warningsDisplayed:
                    print('\nWarning: a reaction of the type "Complex ->  Poly  +  Poly" has been encountered. Although Simply is able to simulate such reactions without any problems, you will have to instruct (through an arm index number) which reactant complex arm becomes which product polymer. By design, the arm index number provided refers to the first polymer product (the other arm will become the second polymer product). Please set up your input accordingly. Note: if scission is possible on more than one place (the case with, or example, RAFT complexes), you will likely need to provide two lines as follows:\n\nComplex[0] -> Poly + RAFT-poly\nComplex[1] -> Poly + RAFT-poly\n\nThis can alternatively be written as:\n\nComplex[0] -> Poly + RAFT-poly\nComplex[0] -> RAFT-poly + Poly\n\n')
                    warningsDisplayed.append(reactionType)
            elif reactionType == 'complex_change':
                if r1_arms != p1_arms:
                    self.reactionTypeError(linecache.getline(file, i), 'Reactant complex and product complex have same number of arms')
                    return -1
            elif reactionType == 'depropagation_complex':
                if r1_arms != p1_arms:
                    self.reactionTypeError(linecache.getline(file, i), 'Reactant complex and product complex have same number of arms')
                    return -1
            elif reactionType == 'complex_cleave_1poly':
                if r1_arms != (p1_arms + 1):
                    self.reactionTypeError(linecache.getline(file, i), 'Product complex should have one arm less than the reactant complex')
                    return -1
                elif r1_arms < 3:
                    self.reactionTypeError(linecache.getline(file, i), 'Reactant complex should have at least three arms')
                    return -1
                elif r1_index >= r1_arms:
                    self.reactionTypeError(linecache.getline(file, i), 'Index of the reactant complex is out of range')
                    return -1
            elif reactionType == 'propagation_complex':
                if r1_arms != p1_arms:
                    self.reactionTypeError(linecache.getline(file, i), 'Reactant complex and product complex have same number of arms')
                    return -1
            elif reactionType == 'complex_expand_1':
                if (r1_arms + 1) != p1_arms:
                    self.reactionTypeError(linecache.getline(file, i), 'Product complex should have one arm more than the reactant complex')
                    return -1
                elif r1_index > r1_arms:
                    self.reactionTypeError(linecache.getline(file, i), 'Index of the reactant complex is out of range')
                    return -1
            elif reactionType == 'complex_combination':
                if (r1_arms + r2_arms) != p1_arms:
                    self.reactionTypeError(linecache.getline(file, i), 'The sum of the number of arms of the two reactants should equal the number of arms of the product')
                    return -1

            # Check if the name is unique
            for j in range(len(moleculesList)):
                if moleculesList[j][1] == line[0]:
                    self.lineError(linecache.getline(file, i),'Names must be unique, however, there is also a molecular species defined with this name')
                    return -1

            for j in range(len(ratesList)):
                if ratesList[j][0] == line[0]:
                    self.lineError(linecache.getline(file, i),'Names must be unique, however, there is also a reaction rate defined with this name')
                    return -1

            # Everything checks out
            reactions.append(line)

        # Check whether there are any molecules that have been specified, but that are not used
        # First, create two lists of molecule names
        mols_in_moleculesList = ['']
        mols_in_reactions = []
        for i in range(len(moleculesList)):
            mols_in_moleculesList.append(moleculesList[i][1])
        for i in range(len(reactions)):
            mols_in_reactions.append(reactions[i][2])
            mols_in_reactions.append(reactions[i][3])
            mols_in_reactions.append(reactions[i][4])
            mols_in_reactions.append(reactions[i][5])

        # For reacting complexes, remove the indices
        for i in range(len(mols_in_reactions)):
            isReactingComplex = self.isReactingComplex('',mols_in_reactions[i],moleculesList)
            if isReactingComplex == -1:
                return -1
            if isReactingComplex[0] == True:
                mols_in_reactions[i] = isReactingComplex[1]
                if mols_in_reactions[i] == -1:
                    return -1

        # Now check the differences between these lists
        mols_difference = set(sorted(mols_in_moleculesList)) - set(sorted(mols_in_reactions))
        if len(mols_difference) != 0:
            print('Warning: the following molecules have been specified by are not used in any reaction: {0}'.format(list(mols_difference)))

        # Check whether there are any rate coefficients that have been specified, but that are not used
        rates_in_ratesList = []
        rates_in_reactions = []
        for i in range(len(ratesList)):
            rates_in_ratesList.append(ratesList[i][0])
        for i in range(len(reactions)):
            rates_in_reactions.append(reactions[i][1])

        # Check differences between these lists
        rates_difference = set(sorted(rates_in_ratesList)) - set(sorted(rates_in_reactions))        
        if len(rates_difference) != 0:
            print('Warning: the following rate coefficients have been specified by are not used in any reaction: {0}'.format(list(rates_difference)))

        # We're done
        return reactions


    def readRateCoeffs(self, file, generalDict, moleculesList):

        # Initialise defaults
        reqLength = [3,4,9,10,6,7,8,10]
        reqTypes = ['f','a','fd','ad','fp','ap','fpr','apr'] # types correspond with above required lengths ('f' has required length 3, etc.)
        rateCoeffs = []

        # Read lines of interest, store each rate in the correct list
        linesList = self.readGroup(file,'rates')
        if linesList == None:
            self.missingGroupError('rates')
            return -1
        elif linesList == -1:
            return -1
        for i in linesList:
            line = linecache.getline(file, i).strip(' \t\n\r').split()

            # Check if the length is OK
            if len(line) not in reqLength:
                self.lineError(linecache.getline(file, i),'Invalid number of arguments')
                return -1

            # Check if the type is known
            line[1] = line[1].lower()
            if line[1] not in reqTypes:
                self.lineError(linecache.getline(file, i),'The rate coefficient is specified incorrectly')
                return -1

            # Check if combination of length and type is OK
            for j in range(len(reqLength)):
                if len(line) == reqLength[j] and line[1] == reqTypes[j]:
                    break
                if j == (len(reqLength) - 1): # Correct combination of length and type is not found
                    self.lineError(linecache.getline(file, i),'Invalid number of arguments')
                    return -1

            # Check for duplicates
            for j in range(len(rateCoeffs)):
                if line[0] == rateCoeffs[j][0]:
                    self.duplicateRateError(line[0])
                    return -1

            # Is the value of the correct type?
            try:
                for j in range(2,len(line)):
                    if (line[1] == 'fpr' and j == 7) or (line[1] == 'apr' and j == 9):
                        pass
                    else:
                        line[j] = float(line[j])
            except ValueError:
                self.valueTypeError('rate coefficient',line[0],line[i],'Could not interpret this value as a floating point number')
                return -1

            # Is the name unique?
            for j in range(len(moleculesList)):
                if moleculesList[j][1] == line[0]:
                    self.lineError(linecache.getline(file, i),'Names must be unique, however, there is also a molecular species defined with this name')
                    return -1

            # Is the temperature required?
            if (line[1] in ['a','ad','ap','apr']) and (generalDict['temperature'] == 0.0):
                if moleculesList[j][1] == line[0]:
                    self.lineError(linecache.getline(file, i),'A temperature dependent rate has been encountered, but no valid temperature has been specified')
                    return -1

            # Should we enable the calculation of the free volume?
            if (line[1] in ['fp','ap','fpr','apr']) and (generalDict['freevolume'] == 0):
                generalDict['freevolume'] = 1

            # Should we enable the recalculation of the conversion at every iteration?
            if (line[1] in ['fd','af','fp','ap','fpr','apr']) and (generalDict['recalcconversion'] == 0):
                generalDict['recalcconversion'] = 1

            # Check if specified molecule names check out
            if line[1] == 'fpr':
                for j in range(len(moleculesList)):
                    if line[7] == moleculesList[j][1]:
                        break
                    if j == (len(moleculesList) - 1):
                        self.lineError(linecache.getline(file, i),'Species {0} is unknown'.format(line[8]))
                        return -1
            if line[1] == 'apr':
                for j in range(len(moleculesList)):
                    if line[9] == moleculesList[j][1]:
                        break
                    if j == (len(moleculesList) - 1):
                        self.lineError(linecache.getline(file, i),'Species {0} is unknown'.format(line[9]))
                        return -1

            # It checks out ok, let's append it
            rateCoeffs.append(line)

        return rateCoeffs


    def readEnthalpies(self, file, generalDict, reactionsList):

        # Initialise defaults
        reqLength = [2]
        enthalpies = []

        # Read lines of interest, store each rate in the correct list
        linesList = self.readGroup(file,'enthalpies')
        if linesList == None:
            if generalDict['simulateheating'] != 0:  # Heating should be simulated, but no enthalpies are specified -> error
                self.missingGroupError('enthalpies')
                return -1
            else:                                    # No enthalpies found, none needed -> return empty list
                return enthalpies
        elif linesList == -1:
            return -1
        if generalDict['simulateheating'] == 0:      # Heating should not be simulated, but enthalpies are specified -> warning
            print('Warning: enthalpies have been specified, but simulation of reaction heating has not been requested')
        for i in linesList:
            line = linecache.getline(file, i).strip(' \t\n\r').split()

            # Check if length of input line is correct
            if len(line) not in reqLength:
                self.lineError(linecache.getline(file, i),'Invalid number of arguments')
                return -1

            # Check if reaction name exists
            for j in range(len(reactionsList)):
                if reactionsList[j][0] == line[0]:
                    break
                if j == (len(reactionsList) - 1):
                    self.lineError(linecache.getline(file, i),'Reaction name not known')
                    return -1

            # Check if the value is correct
            try:
                line[1] = float(line[1])
            except:
                self.lineError(linecache.getline(file, i),'Invalid value specified')
                return -1

            # It checks out ok, let's append it
            enthalpies.append(line)

        return enthalpies


    def readFreeVolume(self, file, generalDict, ratesList):

        # Initialise defaults
        freeVolumeDict = {'vf0':      -1,
                          'tgm':      -1,
                          'tgp':      -1,
                          'alpham':   -1,
                          'alphap':   -1, }
        reqLength = [2]
        reqParams = ['vf0','tgm','tgp','alpham','alphap']

        # Do we need to read this group?
        if generalDict['freevolume'] == 0:
            return freeVolumeDict

        # Read lines of interest, store each rate in the correct list
        linesList = self.readGroup(file,'freevolume')
        if linesList == None:
            self.missingGroupError('freevolume')
            return -1
        elif linesList == -1:
            return -1
        for i in linesList:
            line = linecache.getline(file, i).strip(' \t\n\r').lower().split()

            # Check if length of input line is correct
            if len(line) not in reqLength:
                self.lineError(linecache.getline(file, i),'Invalid number of arguments')
                return -1

            # Check if a parameter by that name is required
            if line[0] not in reqParams:
                self.lineError(linecache.getline(file, i),'Parameter not recognized')
                return -1

            # Check if the parameter was not already specified
            if freeVolumeDict[line[0]] != -1:
                self.lineError(linecache.getline(file, i),'The parameter is specified more than once')
                return -1

            # Check if the value is correct
            try:
                line[1] = float(line[1])
            except:
                self.lineError(linecache.getline(file, i),'Invalid value specified') 
                return -1

            # It checks out ok, let's set the dictionary value and strike the parameter of the requireds list
            freeVolumeDict[line[0]] = line[1]
            reqParams.remove(line[0])

        # Are we done?
        if len(reqParams) != 0:
            self.missingInputError('free volume',reqParams)
            return -1

        return freeVolumeDict


    def checkSpecies(self, line, species, role, moleculesList):
        # If the specified molecule known?

        # Try to look it up directly   
        for j in range(len(moleculesList)):
            if species == moleculesList[j][1]:
                # Check if it should have a index
                if moleculesList[j][0] == 'complex' and role == 'r':
                    self.lineError(line,'The complex "{0}" should have an index number'.format(species))
                    return -1
                else:
                    return True # Found! Let's go back

        # Check if the species is perhaps specified with a index number
        isReactingComplex = self.isReactingComplex(line,species,moleculesList)
        if isReactingComplex == -1:
            return -1
        if isReactingComplex[0] == False:
            self.lineError(line,'Unknown molecule "{0}"'.format(species))
            return -1
        elif isReactingComplex[0] == True and role == 'p':
            print('\nWarning: in the reaction shown below, "{0}" is specified as a product. Since it is a product, the index will be ignored.\n{1}'.format(species, line))

        return True


    def isReactingComplex(self, line, species, moleculesList):

        # Try to understand the species name and stop if we can't
        try:
            # Find opening bracket
            openIndex = [pos for pos, char in enumerate(species) if char == '[']
            if len(openIndex) == 0:
                return (False,species,0)
            elif len(openIndex) == 1:
                pass
            else:
                raise ValueError

            # Find closing bracket
            closeIndex = [pos for pos, char in enumerate(species) if char == ']']
            if len(closeIndex) == 0:
                return (False, species,0)
            elif len(closeIndex) == 1:
                # Complex string should always end with ']'
                if closeIndex[0] != (len(species) - 1):
                    raise ValueError
            else:
                raise ValueError

            # Looks like we've encountered syntax that corresponds with that of a complex
            # Check if it has been defined
            for i in range(len(moleculesList)):
                if moleculesList[i][1] == species[0:openIndex[0]] and moleculesList[i][0] == 'complex':
                    break
                elif i == (len(moleculesList) - 1):
                    self.lineError(line,'Unknown molecule "{0}"'.format(species[0:openIndex[0]]))
                    return -1

            # Is the specified index valid/permitted?
            index = int(species[openIndex[0]+1:closeIndex[0]]) # Raises ValueError in case of failure
            if index > moleculesList[i][2]:
                self.lineError(line,'The index number of "{0}" exceeds its number of arms'.format(species))
                return -1

        except ValueError:
            self.lineError(line,'Error parsing molecule "{0}"'.format(species))
            return -1

        # That's a complex alright. Return a tuple containing the result of the check, the species name, and its index
        return (True,species[0:openIndex[0]],int(species[openIndex[0]+1:closeIndex[0]]))


    def generalizeReaction(self, reaction, moleculesList):
        
        # Initialize
        generalReaction = []
        
        # For each reactant and product in 'reaction', go through 'moleculesList' and look up the molecule type
        # Complexes with indices should be parsed again
        for i in [reaction[2],reaction[3],reaction[4],reaction[5]]:
            isReactingComplex = self.isReactingComplex('',i,moleculesList)
            if isReactingComplex == -1:
                return -1
            if isReactingComplex[0] == True:
                generalReaction.append('complex')
            else:
                for j in range(len(moleculesList)):
                    if i == '':
                        generalReaction.append('')
                        break
                    elif i == moleculesList[j][1]:
                        generalReaction.append(moleculesList[j][0])
                        break

        # Replace each molecule type with its associated value and return the hash
        for n,i in enumerate(generalReaction):
            if i == '':
                generalReaction[n] = 1
            elif i == 'simple':
                generalReaction[n] = 2
            elif i == 'poly':
                generalReaction[n] = 3
            elif i == 'complex':
                generalReaction[n] = 4

        return generalReaction


    def hashReaction(self, reaction):

        # Hash the reaction
        hash = reaction[0]**11 + reaction[1]**7 + reaction[2]**5 + reaction[3]**3
        return hash


    def readGroup(self, file, group):
        
        # Initialize some variables
        no_lines = self.lineCount(file)
        beginFound = False
        endFound = False
        list = []
        
        try:
            # Look through the file; find the first and last line of the group
            for i in range(1,no_lines+1):
                # We're out of lines
                if i == no_lines:
                    return None
                    
                line = linecache.getline(file, i).lower().strip(' \t\n\r')
                # Remove hex characters that are sometimes present at the start of the file
                line = re.sub(r'[^\w#\s]', '', line)
                
                # Are we there yet?
                if beginFound == False and endFound == False:
                    if line == 'begin ' + group:
                        beginFound = True
                elif beginFound == True and endFound == False:
                    if line == 'end ' + group:
                        endFound = True
                        return list
                    elif line[0] == '#':
                        pass
                    else:
                        list.append(i)
                        
        except:
            # All kinds of things could have gone wrong
            print('\nUnknown error while trying to find {0} group'.format(group))
            return -1


    def lineCount(self, file):
        # Count the number of lines in the input file
        with open(file) as f:
            i = sum(1 for _ in f)
        return i+1


    def lineError(self, line, moreInfo):
        # Could not interpret a certain line, so exit
        print('\nError reading the line: "{0}"'.format(line.strip(' \t\n\r').replace('\t',' ')))
        if moreInfo:
            print(moreInfo)
        return -1


    def reactionTypeError(self, line, moreInfo):
        # A certain reaction type is unsupported, so exit
        print('\nError parsing reaction: "{0}"'.format(line.strip(' \t\n\r').replace('\t',' ')))
        if moreInfo:
            print(moreInfo)
        return -1


    def valueTypeError(self, variableType, variable, value, moreInfo):
        # Could not interpret a certain value, so exit
        print('\nError: {0} "{1}" has an invalid value "{2}"'.format(variableType,variable,value))
        if moreInfo:
            print(moreInfo)
        return -1


    def duplicateMolError(self, molecule):
        # Molecule specified more than once, so exit
        print('\nError: molecule "{0}" is specified more than once'.format(molecule))
        return -1


    def missingInputError(self, name, list):
        # Input missing, so exit
        print('\nError: one or more required {0} options have not been specified, namely:\n'.format(name))
        for i in range(len(list)):
            print('   ' + list[i])
        return -1


    def duplicateRateError(self, rate):
        # Molecule specified more than once, so exit
        print('\nError: rate coefficient "{0}" is specified more than once'.format(rate))
        return -1


    def missingGroupError(self, name):
        # Molecule specified more than once, so exit
        print('\nError: Could not find the {0} group. Please check your input file'.format(name))
        return -1