### Simply: a SIMulator for PoLYmerizations

_Current version is 0.98 alpha, released on 2016-06-06 ([download this version](https://github.com/thomaspijper/Simply/releases))_

Simply is a kinetic Monte Carlo simulator for polymerization reactions. It is based on _paraPolySim_ by Barner-Kowollik and coworkers and consists of two separate programs. The first (‘generator’) is Python program that interprets the user-provided input file and translates it into a C header file. The second (‘simulator’) is the actual simulation engine, written in C. The C code, when combined with the generated header file, compiles into an executable that is highly optimized for the system that is to be simulated.

One important goal of the program is to maximize computational efficiency, which is achieved through a) detailed knowledge of the simulation being available at compile time, b) an efficient parallelization scheme (which utilizes MPI), c) a minimalistic code path, and d) support for SSE2 and AVX2 SIMD instructions.

Simply can run on Linux and Windows.

### Capabilities and limitations

The simulation capabilities of Simply are as follows:
* simulation of various types reactions between singles molecules (incl. monomeric species), polymers, and multi-chain polymer complexes;
* simulation of various types of polymerization systems, such as bulk radical polymerization, RAFT, ATRP, etc.
* explicit treatment of the chain lengths of individual polymeric species;
* simulation of heating and cooling of the reaction mixture (by reaction energy or thermal dissipation);
* simulation of autoacceleration kinetics (i.e. the gel effect).

Limitations that should be emphasized are:
* for copolymerizations, Simply does not provide information on the chain sequence of polymers;
* Simply is not designed for on-lattice kMC calculations. For a program that has such capabilities, please have a look at [kMClib](https://github.com/leetmaa/KMCLib).

### Documentation

Documentation is available on the project wiki ([link](https://github.com/thomaspijper/Simply/wiki)).

### Citing
If you’re using Simply for research, please cite it as follows:

    Simply kinetic Monte Carlo simulation package[1], which is based on paraPolySim by 
    Chaffey-Millar et al.[2],[3]

    [1] Simply, by Thomas C. Pijper, https://github.com/thomaspijper/Simply/
    [2] H. Chaffey-Millar, D. Stewart, M. M. T. Chakravarty, G. Keller, C. 
        Barner-Kowollik Macromol. Theory Simul. 2007, 16, 575-592.
    [3] G. Keller, H. Chaffey-Millar, M. M. T. Chakravarty, D. Stewart, C. 
        Barner-Kowollik Specialising Simulator Generators for High-Performance 
        Monte-Carlo Methods Practical Aspects of Declarative Languages 10th 
        International Symposium 2008. P. Hudak, D. S. Warren (Eds.) Springer Berlin 
        Heidelberg, 2008; 4902, 116-132.

### License

Simply’s code is released under the LGPL v3 license. Simply makes use of the following libraries, each of which carries its own license (see the corresponding source code for information):
* [dSFMT](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/) by M. Saito and M. Matsumoto;
* [minunit](https://github.com/siu/minunit) by David Siñuela Pastor;
* [argparse](https://github.com/cofyc/argparse) by Yecheng Fu.

### Contact

Do you have suggestions, bug reports, feature requests, or do you just want to say hello? You can create an _issue_ on the Github page. If this for some reason doesn’t suffice, feel free to drop me a line at tcpijper _at_ gmail _dot_ com. I’m always interested in feedback. As a novice programmer, I’m finding code-related feedback especially valuable.

### Disclaimer

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.
