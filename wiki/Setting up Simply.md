### Linux

Using Simply on Linux systems requires the following software to be installed:

* GCC. Simply has been tested with GCC 5.2, but it will likely work with other versions of GCC as well. In addition, the code should be compatible with the Intel C++ compiler although I'm not able to verify this.
* Python 2.7 (or newer) or 3.5 (or newer). Simply has been tested with Python 2.7.11 and 3.5.1, but older versions of the 2.x and 3.x branches may work as well.
* An MPI implementation. Simply has been tested with OpenMPI, but should work with other MPI implementations as well.

Obtaining and installing these software packages differs which each Linux distribution, so I cannot provide instructions for this.

With the above packages installed, copy the file `Makefile-GCC`in `./simulator/Makefiles/` to `./simulator/` and rename the file to `Makefile`. You can finally make some adjustments to the compiler settings contained in the Makefile, but typically this should not be necessary.

### Windows

Using Simply on Windows systems (Windows 7 or newer) requires the following software to be installed:

* Microsoft Visual Studio 2015. The Community and Express editions ([link](https://www.visualstudio.com/vs-2015-product-editions)) should both work. Which of these versions to use depends on your organization; the use of the Community edition is limited for enterprise organizations. In addition to the Visual Studio C++ Compiler, the code should be compatible with the Intel C++ compiler although I'm not able to verify this.
* Python 2.7 (or newer) or 3.5 (or newer). Simply has been tested with Python 2.7.11 and 3.5.1, but older versions of the 2.x and 3.x branches may work as well. Python can be obtained from [here](https://www.python.org/downloads/).
* An MPI implementation. For Windows, the only MPI implementation I can recommend is Microsoft's. MS-MPI v7.1 can be obtained from the [Microsoft Download Center](https://www.microsoft.com/en-us/download/details.aspx?id=52981). Both packages, _MSMpiSetup.exe_ and _msmpisdk.msi_, should be installed (to their default locations).

With the above installed, copy the file `Makefile-MSVC` in `.\simulator\Makefiles\` to `.\simulator\` and rename the file to `Makefile`. Finally, open the Makefile in a text editor and check if the paths specified at the top of the file are correct. If you installed Visual Studio and MS-MPI to their default folders the paths should be correct, assuming you're using a 64-bit version of Windows.

Finally, go to the Windows _Control Panel_, then _System_. Click on _Advanced System Settings_, then on the _Environment Variables_. Here, edit the `Path` variable and add the line `C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin` (again assuming you installed Visual Studio to its default folder). Note that paths specified here should be separated by a semicolon.

In theory, it should be possible to use MinGW-w64 as a compiler as well. However, setting this up takes some more effort as linking applications to MS-MPI with MinGW-w64 is far from straightforward. I've never attempted to do this myself. I do know that some (outdated?) instructions on how to achieve this float around on the Internet, so if you're interested you may have look at these.