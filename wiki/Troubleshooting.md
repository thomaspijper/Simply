**Simply keeps informing me that the 'monomer audit' has failed**

The purpose of the monomer audit is to verify the consistency of the calculation. This is done by comparing the sum of the lengths of all polymer chains + the number of unreacted monomer particles with the number of monomer particles at the start of the calculation. If Simply finds an discrepancy, it assumes particles got lost during the calculation and informs the user. This may indicate a problem with one of Simply's processes, possibly due to hardware or MPI-related problems. However, if one has specified one or more reactions which consumes monomer particles without increasing polymer chain length, the monomer audit will be incorrect by design. In this case, the message can be safely ignored.

**Simply sometimes reports that the worker time or reduce time is 0 microseconds**

This has been observed on Windows and stems from the fact that on Windows timings are currently retrieved with only millisecond resolution. As such, any time interval shorter than a full millisecond is reported as 0. This will be fixed in the next version of Simply.

**When using MS-MPI, screen output of the different nodes gets mixed up**

This issue has so far only been observed with MS-MPI, not with some of the MPI implementations available on Linux. It is a purely cosmetic issue. A solution is currently being worked on.

**The calculation aborts with the message _warning: state too big (stirr), not yet dealt with_**

This occasionally happens when running calculation in parallel. A solution is being worked on. In the meantime, a workaround for this issue is as follows:

1. Open simply.c
2. Go to line 2345
3. Increase the multiplication factor of _1.2_ to a value where the error message no longer appears. A factor of _1.4_ generally suffices.