### Input structure

Simply's input files structure makes use of groups which contain all data. Each group starts with a line that states `begin <groupname>` and ends with a line that states `end <groupname>`. In between these lines, the rules for input layout are very strict -- not even a blank line is allowed. Outside of groups, the user can write whatever they want as these lines are ignored. Lines of data can contain only alphanumerical characters as well as a few special characters such as the underscore. Words/values may be separated by any number of spaces and/or tabs.

The following group names are accepted by Simply:

* `general` -- used for global options such as the particle count, maximum wall time, switches for certain functionality, etc.
* `molecules` -- used for specifying all molecular species used in the calculation
* `reactions` -- used for specifying all reactions taking place in the calculation
* `rates` -- used for specifying all reaction rate constants
* `freevolume` -- used for specifying parameters used in calculating the free volume fraction of the system (necessary when calculating certain types of rate constants)
* `enthalpies` -- used for specifying reaction enthalpies (for simulating heat generation)

Names of species, reactions, and reaction rates must be unique. A reaction and a rate constant can, for example, not both be named _propagation_.

Values may be given in normal or scientific notation.

An overview of the input accepted for each group is presented below.


### General group

The `general` group accepts only options that are structured as follows:

    <option>   <value>

Available options are as follows:

    particlecount      An integer that specifies the number of particles in the 
                       calculation. This option is required and has no default.

    monomernames       A list of all species in the 'molecules' group that count as 
                       monomer and that is used to differentiate between monomer and 
                       other 'simple' molecules. An addition reaction of a monomer 
                       to a polymer is interpreted as propagation (i.e. chain growth) 
                       while the addition of non-monomeric species (e.g. a RAFT agent) 
                       is not. Multiple species must be separated with commas.

    seed               An integer that is used to seed the random number generator.
                       Changing the seed will cause the result of the simulation to 
                       change. There is no default; in absence of a specified value, a 
                       seed will be chosen at random (the process's PID will be used).

    calcdist           A flag which can be used to control the calculation of the 
                       moments of the distribution, the number-averaged and weight-
                       averaged chain lengths, and polydispersity index. Allowed values
                       are 0 and 1. The default is 0, meaning these parameters are not
                       calculated unless needed for the calculation of rate constants.
                       A value of 1 causes these parameters to be always calculated.

    syncmethod
    syncsimtime
    syncevents         Three keywords that are used to determine the how often 
                       synchronization takes place. During synchronization, particles
                       are redistributed between the different processes (a cyber-
                       equivalent of stirring) and output is written to screen and 
                       files. 'syncmethod' can have one of two values: 'simtime'
                       and 'events'. The default is 'simtime', which causes
                       synchronization to take at a certain interval of simulated
                       time (not wall time). 'events' causes synchronization to
                       take places after a certain number of reaction events have
                       taken place. 'synctime' is an integer that specifies the
                       interval in seconds (default is 5). 'syncevents' is an
                       integer that specified the interval in events (default is
                       500000 or 5e+5).

    maxwalltime        An integer that specifies the maximum computational (wall)
                       time, in minutes. When the wall time has been exceeded,
                       the simulation will be ended after the next synchronization.
                       The default value is 0, meaning the simulation will not end
                       until it is complete.

    maxsimtime
    maxevents
    maxconversion      Three options that specify when the simulation is complete.
                       'maxsimtime' (integer) specifies the maximum simulation 
                       time, 'maxevents' (integer) the maximum number of events,
                       and 'maxconversion' (float) the maximum fraction of
                       monomer that has been used up. The simulation is considered
                       complete when all three criteria have been met. Each has a 
                       default value of 0.

    temperature        A float that specifies the temperature of the simulation, in
                       Kelvin. It also serves as the base temperature to which
                       the system will try to relax by cooling down (or heating up)
                       if a 'coolingrate' is specified. There is no default value.

    coolingrate        The constant 'r' (float) in the formula dT*exp(-rt), used
                       to determine how fast a system relaxes to its base temperature.
                       There is no default value.

    starttemperature   A float that specifies the initial temperature of the 
                       simulation (which may differ from the base temperature).
                       A 'coolingrate' must be specified when using this option.
                       There is no default value.

    heatcapacity       A float that specifies the specific heat capacity of the
                       system, in kJ/(dm^3 K). This is used when simulating reaction 
                       heating.

    longchainsupport   A flag that can be used to increase the maximum allowed
                       chain length from 65,536 to 4,294,967,295. Allowed values
                       are 0 (= 65,536 max) and (= 4,294,967,295 max). The 
                       default is 0. Increasing the maximum allowed chain length
                       will cause Simply to use twice as much system memory.

### Molecules group

Simply considers three types of species, denoted _simple_ for nonpolymeric molecules, _poly_ for polymeric molecules that consist of a single chain, and _complex_ for polymeric molecules that consist of multiple chains. All species that take part in the simulations should be specified in the `molecules` group. Each species should be denoted on a single line, structured as follows:

    simple   <species_name>  <concentration>
    poly     <species_name>
    complex  <species_name>  <number_of_arms>

Hereby, `<species_name>` is a unique name that may consist only of alphanumeric characters and a few special characters (such as the underscore). For _simple_ species, `<concentration>` is an optional float that specifies the initial concentration of the species, in mol/L. When not given, the initial concentration of the corresponding species set as 0. For _poly_ and _complex_ species, the initial concentration is always 0. Finally, `<number_of_arms>` is an integer that specifies the number of polymeric chains that the complex consists of. 

### Reactions group

The `reactions` groups is used to specify all chemical processes. Each line describes a single process and should be of one of the following forms:

    <reaction_name>  <rate_name>  =  <reactant1>                  ->  <product1>
    <reaction_name>  <rate_name>  =  <reactant1>  +  <reactant2>  ->  <product1>
    <reaction_name>  <rate_name>  =  <reactant1>                  ->  <product1>  +  <product2>
    <reaction_name>  <rate_name>  =  <reactant1>  +  <reactant2>  ->  <product1>  +  <product2>

Here, `<reaction_name>` is the unique name of the reaction, `<rate_constant>` is a rate constant defined in the _rates_ group, and reactants and products are species described in the _molecules_  group. The `=`, `+`, and `->` signs should be used verbatim and are used to interpret the line correctly (as well as make it more readable to the user). As noted earlier, any number of spaces or tabs may be used to separate these names and signs. A line may not contain more than two reactants and two products.

Whenever a _complex_ species is used as a reactant, an arm index should be specified by adding it to the end of the species name, in between square brackets. An example:

    begin molecules
      poly     P   
      poly     P_RAFT
      complex  Q        2
    end molecules
    
    begin reactions
      fragmentation1   k_frag   =   Q[0]     ->    P    +    P_RAFT
      fragmentation2   k_frag   =   Q[1]     ->    P    +    P_RAFT
    end reactions

The arm index has the following meanings:

    Complex  +  None     ->  Poly     +  None
    # No function
    
    Complex  +  None     ->  Poly     +  Poly
    # Arm index specifies the chain to become the first Poly product
    
    Complex  +  None     ->  Complex  +  None
    # No function
    
    Complex  +  None     ->  Complex  +  Simple
    # Arm index specifies the chain from which the Simple species is cleaved
    
    Complex  +  None     ->  Complex  +  Poly
    # Arm index specifies the chain which dissociates from the complex
    
    Complex  +  Simple   ->  Complex  +  None
    # Arm index specifies the chain to which the Simple species is added
    
    Complex  +  Poly     ->  Poly     +  None
    # No function
    
    Complex  +  Poly     ->  Complex  +  None
    # Arm index specifies the position at which the chain in inserted
    
    Complex  +  Complex  ->  Poly     +  None
    # No function
    
    Complex  +  Complex  ->  Complex  +  None
    # No function

Note: because of how input files are parsed, the arm index should _always_ be specified, even if it has no purpose. This may be improved upon in a future release.

### Rates group

The `rates` group is used to specify the rate constants of the reactions specified in the _reactions_ group. Depending on what kind of rate constant is specified (see the part of the wiki on Theory), input should be of one of the following forms:

    <rate_name>  f    <value>
    <rate_name>  a    <A>  <Ea>
    <rate_name>  fps  <value>  <kd,0>  <A>  <B>
    <rate_name>  aps  <A>  <Ea>  <kd,0>  <A>  <B>
    <rate_name>  fd   <value>  <d1>  <d2>  <d3>  <d4>  <d5>  <d6>
    <rate_name>  ad   <A>  <Ea>  <d1>  <d2>  <d3>  <d4>  <d5>  <d6>

Here, the second item in the line specifies what type of rate constant is specified. The following types are available:

* Type `k` specifies a fixed value. `<value>` should be in either 1/s or L / (mol s), depending on the reaction order.
* Type `a` specifies a rate constant with Arrhenius parameters. `<A>` is the pre-exponential factor and should be in either 1/s or L / (mol s), depending on the reaction order. <Ea> is the activation energy, in kJ/mol. The rate constant is calculated as follows:
    
    ![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula11.png)

* Type `fps` specifies a rate constant that consists of a 'chemical' part and a 'diffusion' part according to the formula:
    
    ![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula12.png)
    
    with the 'diffusion' part being calculated as follows:
    
    ![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula13.png)
    
    `<value>` is a fixed value for _kc_ that should be in either 1/s or L / (mol s), depending on the reaction order. `<kd,0>` is value fo _kd,0_, also in either 1/s or L / (mol s). `<A>` and `<B>` are the values for parameters _A_ and _B_. The use of this type of rate constant also requires the user to specify parameters used to calculate the fractional free volume (see the _freevolume_  group).

* Type `aps` is identical to type `fps`, except for the 'chemical' part of the rate constant (_kc_) which is defined in terms of Arrhenius parameters.

* Type `fd` consists of a fixed value for _kc_ multiplied by a six-parameter exponential decay function:
    
    ![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula15.png)
    
    Hereby, values of `<d*>` correspond to the six parameters used in this function.

* Type `ad` is identical to type `fd`, except for the rate constant being defined in terms of Arrhenius parameters.

### Freevolume group

The `freevolume` group is used to specify parameters used for calculating the fractional free volume. The formula used for this calculation is as follows:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula14.png)

with _X_ being the conversion.

Input for this group should be of the following form:

    <parameters>   <value>

Available parameters are:

    vf0        The constant vf,0 in the above formula, which  is 
               is the fractional free volume when the temperature
               of the system is equal to its glass transition 
               temperature.
    tgm        The glass transition temperature of the system at
               t = 0, in K.
    tgp        The glass transition temperature of the system at
               complete conversion (i.e. the polymer), in K.
    alpham     The volume expansion coefficient of the system at
               t = 0, in 1/K.
    alphap     The volume expansion coefficient of the system at
               complete conversion, in 1/K.

All of these parameters should be specified as none of them has a default value.

### Enthalpies group

The `enthalpies` group is used to define enthalpies of reaction for reactions specified in the _reactions_ group. Input for this group should be of the following form:

    <reaction_name>  <value>

where `<reaction_name>` is the name of the reaction as defined in the _reactions_ group and `<value>` is the enthalpy of reaction for this reaction, in kJ/mol. If a specified reaction is not named in this group, its enthalpy of reaction is assumed to be 0.