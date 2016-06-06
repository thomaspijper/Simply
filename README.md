### Simply: a SIMulator for PoLYmerizations

Simply is a kinetic Monte Carlo simulator for polymerization reactions. It is based on _paraPolySim_ by Barner-Kowollik and coworkers and consists of two separate programs. The first (‘generator’) is Python program that interprets the user-provided input file and translates it into a C header file. The second (‘simulator’) is the actual simulation engine, written in C. The C code, when combined with the generated header file, compiles into a executable that is highly optimized for the system that is to be simulated.

One important goal of the program is to maximize computational efficiency, which is achieved through a) detailed knowledge of the simulation already being available at compile time, b) an efficient parallelization scheme (which utilizes MPI), and c) a code design that minimizes the amount of operations needed for each iteration.

Simply can run on Linux and Windows.

### Capabilities and limitations

The simulation capabilities of Simply are as follows:
* simulation of various types reactions between singles molecules (incl. monomeric species), polymers, and multi-chain polymer complexes;
* simulation of various types of polymerization systems, such as bulk radical polymerization, RAFT, ATRP, etc.
* explicit treatment of the chain lengths of individual polymeric species;
* simulation of heating (by enthalpy of reaction) and cooling of the reaction mixture;
* the ability for the user to specify rate constants that depend on additional variables, such as temperature, free volume, conversion, etc.

Limitations that should be emphasized are:
* for copolymerizations, Simply does not provide information on the chain sequence of polymers;
* Simply cannot perform on-lattice kMC calculations. For a program that has such capabilities, please have a look at [kMClib](https://github.com/leetmaa/KMCLib).

### Documentation

Documentation is available on the project wiki ([link](https://github.com/thomaspijper/Simply/wiki)).

### Citing
If you’re using Simply for research, please cite it as follows:

    Simply Monte Carlo simulation package[1], which is based on paraPolySim by 
    Chaffey-Millar et al.[2],[3]

    [1] Simply, https://github.com/thomaspijper/Simply/
    [2] H. Chaffey-Millar, D. Stewart, M. M. T. Chakravarty, G. Keller, C. 
        Barner-Kowollik Macromol. Theory Simul. 2007, 16, 575-592.
    [3] G. Keller, H. Chaffey-Millar, M. M. T. Chakravarty, D. Stewart, C. 
        Barner-Kowollik Specialising Simulator Generators for High-Performance 
        Monte-Carlo Methods Practical Aspects of Declarative Languages 10th 
        International Symposium 2008. P. Hudak, D. S. Warren (Eds.) Springer Berlin 
        Heidelberg, 2008; 4902, 116-132.

### License

All code is released under the LGPL v3 license, with the exception of the dSFMT pseudorandom number generator (written by M. Saito and M. Matsumoto) which is released under the 3-clause BSD license.

### Contact

Do you have suggestions, bug reports, feature requests, or do you just want to say hello? You can create an _issue_ on the Github page. If this for some reason doesn't suffice, feel free to drop me a line at tcpijper _at_ gmail _dot_ com. I’m always interested in feedback. As a novice programmer, I'm finding code-related feedback especially valuable.

### Disclaimer

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.
