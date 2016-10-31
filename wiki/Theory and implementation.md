### Simply’s kinetic Monte Carlo method

Simply starts with an analysis of all specified species, reactions, and rates constants. Based on the concentrations of species and specified rate constants, for each reaction the chemical rate is calculated. For a unimolecular reaction, the rate is simply:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula1.png)

with _XA_ being the number of particles of a species _A_. For a bimolecular reaction, the chemical rate equals either:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula2.png)

in the case of non-identical reacting species, or:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula3.png)

in the case of identical reacting species. As the rate constant _ki_ must be in units of _1/s_, it is modified in the case of a bimolecular reaction according to the following formula:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula4.png)

whereby an additional multiplication by 2 is required in the case of identical reacting species. It should be noted that these transformations are automatically made by Simply’s generator program and should not be performed by the user.

Each chemical rate _Ri_ is then converted to a probability _Pi_ according to the formula:

   ![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula5.png)

with _r_ being the total number of reactions.

From this initial state, the Monte Carlo algorithm proceeds as follows:
<br/>
<br/>

1) A random number _u_ is calculated by multiplying a random number on the interval _(0,1]_ by the cumulative probability ![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula6.png). For a reaction with an index _n_ in the range _1 = n = r_, if _u_ lies in the range:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula7.png)

reaction _n_ is the reaction to be performed.

2) The chosen reaction is performed and the particle counts of the species involved are adjusted. If the reaction involves one or two polymeric species, the reacting chain lengths are chosen using new random numbers on the interval _(0,1)_. (Note: formally, an interval of [0,1] should be used, but this is not possible due to limitations of the PRNG. It is expected that an interval of (0,1) produces virtually identical results.)

3) The time _t_ is advanced according to the formula:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula8.png)

with _u’_ being a new random number on the interval _(0,1]_.

4) The initial state is recalculated in order to account for change in particle counts. Simply will hereby only recalculate the rates of reactions for which the particle counts of its reactant(s) were changed. The algorithm will then proceed with step 1.
<br/>
<br/>

In an effort to perform step 1 as efficiently as possible, Simply does not perform a linear scan through the different reaction probabilities until a reaction _n_ has been identified. Instead, the reaction is chosen through a binary probability tree with a number of leaves that a) is a power of 2 and b) that is at least equal to _r_. As an example, for a system consisting of 5 reactions, the probability tree is as follows:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/tree.png)

Excess leafs hereby carry a probability of 0. The use of a binary tree gives the selection of a reaction a time complexity of only _O(log r)_ while a linear scan would have a time complexity of _O(r)_.

### Supported reaction types

Simply distinguishes three types of species, denoted _simple_, _poly_, and _complex_. Hereby, a _simple_ species is a non-polymeric molecule, a _poly_ species is a linear polymer chain, and a _complex_ species is a polymer consisting of more than one chain. With these three types, the following reactions are possible:

    Simple               ->  Simple
    Simple               ->  Simple   +  Simple
    Simple   +  Simple   ->  Simple
    Simple   +  Simple   ->  Simple   +  Simple
    Simple               ->  Poly
    Simple   +  Simple   ->  Poly
    Poly                 ->  Poly
    Poly                 ->  Poly     +  Simple
    Poly     +  Simple   ->  Poly
    Poly     +  Simple   ->  Poly     +  Simple
    Poly     +  Simple   ->  Poly     +  Poly
    Poly     +  Poly     ->  Poly
    Poly     +  Poly     ->  Poly     +  Poly
    Poly     +  Poly     ->  Complex
    Complex              ->  Poly
    Complex              ->  Poly     +  Poly
    Complex              ->  Complex
    Complex              ->  Complex  +  Simple   
    Complex              ->  Complex  +  Poly
    Complex  +  Simple   ->  Complex
    Complex  +  Poly     ->  Poly
    Complex  +  Poly     ->  Complex
    Complex  +  Complex  ->  Poly
    Complex  +  Complex  ->  Complex

### Rate constants

Various types of rate constants can be specified by the user. The most straightforward is to specify a single, constant value or to specify the rate constant in terms of the Arrhenius parameters _Ea_ (the activation energy, in kJ/mol) and A (the pre-exponential factor, in either 1/s or L/(mol s) depending on the reaction). In the latter case, the rate constant is calculated according to the Arrhenius equation:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula11.png)

Here, _R_ is the gas constant (in J / (mol K)) and T is the temperature (in K).

It is also possible to specify rate constants that take into account the autoacceleration effect (commonly referred to as the gel effect or Trommsdorff effect) that is observed with bulk polymerizations. One such rate constant implemented in Simply is of the form:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula12.png)

where _kc_ and _kd_ are the 'chemical' and 'diffusion' contributions to the rate constant _k_, respectively. The former is hereby either a fixed value or calculated from Arrhenius parameters, as described above. The latter is calculated as follows:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula13.png)

where _kd,0_ is a constant value (in either 1/s or L/(mol s) depending on the reaction), _Xn_ is the number average chain length, _vf_ is the fractional free volume, and _A_ and _B_ are additional parameters. The fractional free volume _vf_ is calculated as follows:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula14.png)

Here, _vf,0_ is a constant, _alpha_ is the volume expansion coefficient (in 1/K), _T_ is the temperature (in K), _Tg_ is the glass transition temperature (in K), and _x_ is the conversion (which ranges from 0 to 1). Subscripts _M_ and _P_ refer to monomer and polymer, respectively.

It is also possible to simulate the autoacceleration effect is through to the following formula:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula15.png)

which calculates the rate constant based on six empirically derived parameters as well as the conversion _X_ and (optionally) the temperature _T_. Here, _kc_ is the 'chemical' rate constant which is either a fixed value or calculated from Arrhenius parameters, as described above.

Note: rate constants should always have a unit of 1/s for unimolecular reactions, and L/(mol s) for bimolecular reactions.

### Heating and cooling

Simply is able to simulate heating of the system by reaction heat. For this, the user should provide the enthalpy of reaction for any of the reactions specified, as well as the specific heat capacity of the system. Heating is simulated by using the following formula:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula9.png)

where _deltaHr_ is the enthalpy of reaction (in kJ/mol), _C_ is the specific heat capacity (in J/(K L)), _NA_ is the Avogadro constant, and _V_ is the volume of the system (in L).

Cooling of the system can also be simulated. For this, the following exponential decay formula is used:

![](https://github.com/thomaspijper/Simply/blob/master/wiki/formula10.png)

Here, the constant _r_ (the cooling rate in 1/s) should be specified by the user. The time step used is the one calculated above in step 3) of the Monte Carlo algorithm.

Note: it is not possible for the user to define a temperature gradient.

### Parallelization of simulations

Parallelization is achieved by splitting up the simulated system in multiple smaller systems, which are synchronized ('mixed') at user-defined intervals. A simulation is structured as follows:

* When the simulation is initialized, the particles of each species are equally divided among the participating nodes. Reaction rates are adjusted accordingly to account for the reaction volume being divided up.
* Each node runs individually until the user-defined synchronization interval has elapsed.
* At synchronization, each node packages its data (particle counts, chain lengths, simulation time, simulation temperature, etc.).
* The package is communicated to all other nodes, thus providing each node with the full simulation state.
* All nodes take a part of each species particle count. Simulation time and temperature are averaged.
* All nodes recalculate their probability tree.
* Each node continues it simulation until the next synchronization point.

Before and after synchronization, the amount of reacted monomer is compared with the sum of all chain lengths to ensure the consistency of the simulation.

Because each node will simulate a small part of the system, with regular 'mixing' between nodes, a simulation executed on a single node will provide different results than the same simulation run on multiple nodes. For simulations executed on a large number of nodes, regular synchronization is required to obtain more accurate results, i.e. results that are closer to those that would be obtained with a single node. More accurate results could of course also be obtained by increasing the total particle count of the system, but the accompanying increase in computation cost would then defeat the purpose of parallelization. As such, in order to obtain accurate results, a balance between particle count, number of nodes, and synchronization interval will have to be found. The [paper by Barner-Kowollik and co-workers](https://www.cse.unsw.edu.au/~chak/project/polysim/) provides some good pointers on this. Note: on SMP systems, the synchronization step is expected to contribute significantly less to the calculation time compared to systems on which nodes do not share their memory.