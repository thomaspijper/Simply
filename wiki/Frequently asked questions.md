**How can I choose a suitable system size (particle count)?**

That is probably the most often asked question in kinetic Monte Carlo simulations. My advice is to observe the least occurring reactions and least numerous species. If these are not significantly affected by an increase in system size, you _likely_ have chosen a size large enough to lead to accurate results.

**What is the optimal synchronization interval?**

First, it is important to restate that synchronization is also the time at which output is written to the output files. So, if you would like to observe the evolution of your system in detail, you should set the synchronization interval small.

For nodes which make use of the same system memory (SMP systems), I expect the overhead of the synchronization step to be rather small, so a small sync interval should not impair computational speed significantly here. Frequent synchronization will make results more accurate, so for these systems it might be better to sync often. However, for nodes which do not share memory (such as those that communicate over InfiniBand), I expect the overhead to be larger, especially with large particle counts. For such systems, excessive synchronization may have a detrimental on computational speed. I advice you to perform some tests with your particular system. In addition, the [paper by Barner-Kowollik and co-workers](https://www.cse.unsw.edu.au/~chak/project/polysim/) (which describes the design of paraPolySim, on which this project is based) provides some good pointers on this. 

**How can I simulate reactions involving more than two reactants or products?**

Simulating a reaction that gives three products is easily done by splitting the reaction up in two reactions and using an intermediate species with an extremely short lifetime. 

    <reactant1>  +  <reactant2>  ->  <product1>  +  <product2>  +  <product3>

can be written as:

    <reactant1>  +  <reactant2>  ->  <product1>  +  <productX>
    <productX>                   ->  <product2>  +  <product3>

whereby the second reaction has a very high reaction rate (causing the intermediary _productX_ to be used up immediately after formation).

A reaction with three reactants takes more effort though. In reality, trimolecular reaction are very rare, so you should be able to write up any such reaction as two separate bimolecular reactions. However, adjusting the rate constants accordingly is, I believe, complicated.

**Is Simply capable of simulating other chemical processes?**

Simply does not exclusively have to be used for simulating polymerizations. However, you will always need to specify one or more species as monomer. This does not affect results though as this information is only used for the monomer audit. During simulation, the monomer audit will likely give false feedback, but this can be safely ignored.

**Does Simply benefit from Hyper-threading?**

It is reported that Simply benefits from Hyper-threading, though the exact speed-up depends on your processor type.

**Can I compile Simply as a 32-bit application instead of a 64-bit one?**

Although I do not support it, compiling the simulator program as 32-bit will likely work. However, the explicit storage of chain lengths typically requires a large amount of memory, so you may run into memory limitation problems quite quickly.

**Will Simply run on UNIX-based systems other than Linux?**

I expect Simply to run on BSD and OS X/macOS systems as well but have never tried. Please let me know if you have tested this.