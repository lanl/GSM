# The Generalized Spallation Model and Event Generator
___
The Generalized Spallation Model, herein GSM, is a *mostly* object-oriented event generator for spallation events.  
`Physics`: GSM is composed of several sub-models, consisting primarily of the Standard and Modified Dubna Cascade Models (DCM), a Coalescence model, Preequilibrium model, and the GEM2 (Generalized Evaporation Model), being used for simulating compound Evaporation and fission.  
`Code`:GSM utilizes primarily Fortran2003 and Fortran2008, however does have some Fortran66 and Fortran77 within its Modified DCM sub-model.  
`Parallel`: The GSM event generator is parallelizable for single-event simulations when **not** using the `Modified DCM` sub-model at this time.  

> **GSM Version 1.0.0 beta**


## Table of Contents
___
1. [Disclaimer](#disclaimer)
2. [License](#license)
3. [Build Requirements](#build-requirements)
4. [Compilation](#compilation)
   -  [Testing](#testing)
5. [Usage](#usage)
   - [Standalone Usage](#standalone-usage)
   - [Client Usage](#client-usage)
6. [Documentation](#docs)
7. [Contribute](#contribute)
8. [Credits](#credits)


## Disclaimer <a name="disclaimer"></a>
___
>NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL SECURITY LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.



## License <a name="license"></a>
___
The predecessor of GSM, CEM, has been approved for release with associated LA-CC number LA-CC-04-085. The GSM event generator has **not** been licensed at this time for public release.
The CEM material was procuded under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory for the U.S. Department of Energy. The Government is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this material to reproduce, prepare derivative works, and perform publicly and display publicly.  

**The GSM event generator code has not, at this time, been approved for public release under the CEM license.**  


## Build Requirements <a name="build-requirements"></a>
___
GSM was migrated to an object-oriented framework using primarily GNU Fortran with GCC-6.3.0. The following table shows the compilers that can successfully compile GSM:

| Compiler         | Version       | Note                                                  |
| :--------------- | :------------ | :---------------------------------------------------- |
| GCC              | 4.6.1         | Regression differences seen with GCC4.7.2 and earlier |
| Intel            | Not Tested    | |
| PGI              | Not Tested    | |
| LLVM             | Not Tested    | |
| Cygwin           | Not Tested    | |


## Compilation <a name="compilation"></a>
___
A GSM executable can be created by navigating to the primary GSM top-level directory with any of the compilers supported (see the [build requirements](#build-requirements) section for more information). GSM uses a relatively basic CMake build system, **requiring** out-of-source builds. To install GSM, do the following:
```bash
mkdir <build-src>
cd <build-src>
cmake <options> <path-to-GSM-top-level-directory>
make
```

More explicitly stated, users type the following to build GSM in the top-level directory of GSM:  
```bash
mkdir ./build
cd ./build
cmake ../
make
```

> Note: the CMake build platform was created using Linux; this CMake distribution has NOT been tested on non-linux operating systems nor has it been tested.


The GSM installation is done simply by creating a directory outside the top-level directory and `./src/` subdirectories of GSM. In the above example for the out-of-source build, users are recommended to create a directory called `build`, say in the top-level directory of GSM. Users then `cd` to the directory and run CMake (``cmake ..``), generating all of the required Makefiles for each of the GSM libraries. To then compile GSM, users can type ``make`` in the `./build/` directory to generate all necessary object, module, and sub-module files.  
All compiled module and sub-module files are stored in the `build` directory's `./modules/` folder.  
All libraries created by `Make` during compilation are stored in the `build` directory's `./lib/` folder.  
The executible for GSM, `xgsm1`, is housed in the `build` directory's `./bin/` folder.

> GSM requires several data files and tables for its simulations to be in the directory where simulations are done. Upon being built, GSM provides such a directory, called `my_GSM`, in the top-level directory of GSM. This folder contains all required data tables and files in addition to the most recently built version of GSM. Note that this folder is *only* created when GSM is built successfully.


### Testing <a name="testing"></a>
Upon creation of the GSM executable, users should perform all provided regression tests to verify the proper compilation f the GSM event generator. The folder `test`, housed in the top-level GSM directory, contains all available tests for users.  
A script is provided to perform all regression tests to compare the results produced by GSM when built by `GCC6.3.0`. A limited number of tests are currently provided (30 at this time) with small event limits. It is recognized that the provided regression suite is not particularly thorough, however it is assumed to be sufficient.  
Users can test their installation of GSM using the provided regression script from the top-level directory of GSM by typing the following:
```bash
sh -e test/regression/bin/regression.bc
```

The script will attempt to create a clean installation of GSM by creating a `./build/` directory, compiling GSM, and performing all regression tests.

> NOTE: The regression script utilizes a Python module (called `Comparer.py`) to compare output files from the previous the current versions of GSM. This python script will create a file for *each* set of input files tested stating the differences in the output file and flagging them. This Python script was developed utilizing `Python3.6`. It is possible that regression tests will successfully complete with an earlier version of Python, however it is not guaranteed to work properly.

Regression tests can be updated easily with the following command, again being typed in the top level directory of GSM:
```bash
sh -e test/regression/bin/updateOutputs.bc
```

This script will simply copy all files from the `./resultsCurrent/` directory to the `./resultsPrevious/` directory and change the files' names. This script was found helpful when either (a) increasing the number of events simulated for an input file, (b) when adding input files for testing, or (c) during the development of GSM for discovered bugs, modifications, improvements, or any combination of these.



## Usage <a name="usage"></a>
___
Upon being compiled and testing (see the [Compilation](#compilation) and [testing](#testing) sections), users, clients, or both can create and utilize GSM for their spallation simulations.


### Standalone (User) Usage <a name="standalone-usage"></a>
GSM, as a standalone code, is contained by the `standaloneGSM` [module](src/DriverGSM/standaloneGSM.f90). An exectuble, named `xgsm?`, where the `?` is the current major version of GSM, is created upon successful [compilation](#compilation) in a directory named `my_GSM`, named accordingly in the [CMakeLists.txt file](./CMakeLists.txt).

> Users start a simulation by typing `./xgsm?`, where the `?` is replaced by the current MAJOR version number according to the [CMakeLists.txt file](./CMakeLists.txt).

Recent improvements to GSM allow users to specify command line arguments. The following table details the available arguments:  

| Command    | Argument                            | Default Value      | Details                                                               |
| :--------- | :---------------------------------- | :----------------- | :-------------------------------------------------------------------- |
| -i=        | Input file Name                     | `gsm1.inp`         | The name of the input file, including or excluding path name. If spaces are present, enclose in quotation marks |
| -v=        | Verbosity Number (integer or real)  | N/A                | Allows end-user to set the simulation verbosity                       |
| -help      | N/A                                 | N/A                | Prints all valid command line arguments to the end-user               |

Prior to simulation, the GSM driver verifies that all arguments given are valid and checks for file existence(s). Arguments can be given to the GSM version 1 executable, for example, as `./xgsm1 -v=1 -i=myInputFile.inp`. The driver will then read the input file, given the file exists, and proceed to perform all simulations specified within the input file. Information about the simulation will be printed to the user in the form of comments, warnings, errors, and in rare cases, fatal errors, depending on the verbosity set by default or by the user.

> Recommendation: A verbose filter of two or three for standalone simulations for users, and up to four or five for developers, is recommended.
  
> Tip: Port all output to another file during the simulation by appending `> {someFileName}`, with `{someFileName}` being the desired file name, to the command line. Doing so will port all information printed to `standard_out` to the file. Information printed to the `standard_error` unit will still appear on the terminal during the simulation. Information printed to the terminal **will** be more useful for diagnosing problems encoutered when using GSM, particularly for developers.


Given a valid input file name, the GSM driver will proceed to read the input file and perform all simulations specified within the input file.
The `CEM03.03` Manual, [LA-UR-12-01364](https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-ur-12-01364.pdf), provides a fairly in-depth description of the `CEM03.03` input file, being fairly similar to that of GSM. The GSM input file is similar, however with the following modifications and additions:
1. A ZAID (ZZZAAA) may be used to specify light- and heavy-ion interactions. In instances where a non-zero quantum decay number is desired, users enter the incident particle as `Z, A, S`, where `Z` and `A` are the atomic and mass numbers, respectively, and `S` is the desired quantum decay number.
2. The `npreqtyp` card specifying the number of particles considered for preequilibrium emission is inserted immediately proceeding the `nevtype` card.
3. A `pisa` card was added immediately proceeding the `npreqtyp` card. The card specifies whether or not to calculate double differential [mb/sr/MeV], angle integrated [mb/MeV], and energy integrated [mb/sr] spectra for all particles considered, ranging from `n` to `Mg-28`.
   - If the `pisa` card is greater than zero, the `pisa` spectra described will be printed in the output file. Users *must* provide, on the proceeding line, ten base angles for the calculation in addition to a bin width, being constant for each bin. The bin width is recommended to be `5.0`.
4. The `ihist` card was added proceeding the `pisa` card, and the appropriate angles if used, to flag whether or not to print histogram information on the residual nuclei created during the simulation.


### Client Usage (GSM API) <a name="client-usage"></a>
What is up with the API? Talk about data initialization, model construction, etc.


## Documentation <a name="docs"></a>
___
Documentation of the GSM prior to its creation may be found by referencing various CEM and LAQGSM event generator reports. Note that some of the documentation contained by GSM refers to some of these reports and the references therein. Documentation of the source code for GSM may be easily generated via Doxygen as follows:

```bash
cd ./Docs
doxygen Doxyfile.in
```

Doing so generates, by default, documentation in both a `LaTeX` and `HTML` format. The generated `HTML` documentation may be viewed by pointing a browser to the `index.html` file within the `Docs/Docs_HTML` directory.


## Contribute <a name="contribute"></a>
___
Contributions to GSM can be made by contacting Dr. Leslie Kerby (kerblesl@isu.edu) at the Idaho State University or by contacting the [XCP-3](https://www.lanl.gov/org/padwp/adx/computational-physics/monte-carlo/index.php) group leader, currently Avneet Sood (sooda@lanl.gov), of the Los Alamos National Laboratory with information regarding the contribution.  

>Tip: developers can obtain a general understanding of the code styling within GSM by looking through the [Coalescence model code](./src/Coalescence). This code represented there is sufficiently short and clear to provide beginning developers with style reference and the architecture utilized for the GSM code.  

Questions about GSM, regarding code or physics, can be directed to C. Juneau (junechas@isu.edu) and Dr. Leslie Kerby (kerblesl@isu.edu).


## Credits <a name="credits"></a>
___
The primary authors of GSM, originally CEM, include Drs. S. Mashnik, K. Gudima, and A. Sierk, with important contributions from R. Prael, M. Baznat, and N. Mokhov. Upgrades made to CEM03.03 prior to GSM's creation can be mostly accredited to Dr. Leslie Kerby, who provided several important updates to the preequilibrium model of GSM. GSM was migrated to an object-oriented framework by C. Juneau under the auspices of the Idaho State University and, in part, from funding provided by the Los Alamos National Laboritory under Idaho State University subcontract number 385443.
