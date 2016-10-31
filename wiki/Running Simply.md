### Using the generator program

The generator program converts the plain text input file into a C header file. Using it is simple: just run the generator with the `generateCode.py` file and input the name of the input file. If the input file does not reside in the same folder as `generateCode.py`, the path should be specified as well. The generator program will then create a file named `genpolymer.h`, which is the C header file required to compile the simulator program.

`generateCode.py` accepts the following optional command line options:

    -i <input_file>      Specifies the name (optionally including the path)
                         of the input file.
    -i <output_file>     Specifies the name (optionally including the path)
                         of the output file. Note that the file should be 
                         ultimately be called genpolymer.h, otherwise the
                         simulator program will not compile successfully.
    -f                   Forces the output file to be overwritten if it
                         already exists.
    -v                   Verbose mode. The generator program will print
                         its interpretation of the input file.
    -h                   Prints information about the available command
                         line options.

### Compiling the simulator program

**Linux**

On Linux systems, the simulator program can be compiled using the command `make` (which should be given while in the `./simulator` folder). This results in the creation of an executable named `simply`. Once the simulation is complete and the output files have been copied elsewhere, the command `make clean` can be used to delete the binary and the `.o` object files that were created during compilation. Important: this also deletes the `.csv` output files!

**Windows**

Compilation on Windows requires a little bit more effort compared to compilation on Linux. You can compile using Visual Studio's developer command prompt, however, I personally prefer to compile from the normal command prompt. Open a command prompt, navigate to `C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC`, then issue the command `vcvarsall.bat amd64`. This sets all the environmental variables required for compiling on the command line. The command will have to issued each time you (re)open the command prompt.

From here on, compiling is straightforward. While in the `simulator` folder, type `nmake` to compile the simulator program. Once the simulation is complete and the output files have been copied elsewhere, the command `nmake clean` can be used to delete the binary and as well as various intermediary files that were created during compilation. Important: this also deletes the `.csv` output files!

### Using the simulator program

Running the simulator program as a single process is achieved by directly running the executable. Running the simulator program in parallel should however be done through the `mpiexec` or `mpirun` program (depending on which MPI implementation you installed). A typical command is as follows:

    mpiexec -np 4 simply > output.txt 2>&1

This starts a run on 4 nodes and writes all screen output to the file `output.txt`. The addition of `2>&1` causes error messages to be written to the output file as well.

The simulator program accepts optional command line options that can be used to adjust some parameters of the simulator program (without the need to recompile). These are as follows:

    -s <seed>            Specifies a seed for the pseudorandom number
                         generator. <seed> should be an integer.
    -e <events>          Specifies the synchronization interval in terms
                         of reaction events. This option only has an 
                         effect if 'syncmethod' was set to 'events'. 
                         <events> should be an integer.
    -s <simtime>         Specifies the synchronization interval in terms
                         of simulation time. This option only has an 
                         effect if 'syncmethod' was set to 'simtime'. 
                         <simtime> should be a float.
    -h                   Prints information about the available command
                         line options.

_Note: command line options are currently only supported on Linux._


### Output

Output files written by the simulator program are in the csv format. Files that are written during simulation and that can be used to follow the evolution of the system are `concentrations.csv` and `rates.csv`. The former file contains the concentrations of all species as well as the conversion. Optionally, the temperature and moments of the distributions (including averaged chain lengths and polydispersity index) are written as well. The latter file contains the rate constants of all reactions.

Upon completion of the simulation, Simply writes detailed information regarding particle counts and chain lengths for each polymeric species. The corresponding files are named `<name>-<time>.csv` where `<name>` is the name of the species and `<time>` is the simulation time at which the simulation reached completion.