# The Generalized Spallation Model and Event Generator

The Generalized Spallation Model, herein GSM, is a *mostly* object-oriented event generator for spallation events. GSM is intended to simulate particle interactions for energies greater than _100 MeV/A_, up to several _TeV/A_.

> **GSM version 1.0.0 beta**

## Contents

+ [License and Disclaimer](#license-and-disclaimer)
+ [The GSM Event Generator](#about-the-gsm-event-generator)
   + [Simulation Flow](#simulation-flow)
   + [Software Architecture](#software-architecture)
+ [Build Requirements](#build-requirements)
+ [Compilation](#compilation)
   + [Testing](#testing)
+ [Usage](#usage)
   + [Standalone Usage](#standalone-usage)
   + [Client Usage](#client-usage)
+ [Documentation](#docs)
+ [Contribute](#contribute)
+ [Credits](#credits)



## License and Disclaimer <a name="license-and-disclaimer"></a>

The GSM is an open source program under the BSD-3 license. See the [LICENSE](LICENSE) and [COPYRIGHT](COPYRIGHT) for more details.


## About the GSM Event Generator <a name="about-the-gsm-event-generator"></a>

The Generalized Spallation Model is a Monte-Carlo event generator meant to simluate particle-particle collisions ranging from ~_100 MeV/A_ up to several _TeV/A_. The GSM utilizes a modern software architecture to help facilitate parallelized model physics simulations, as well as providing interchangeable interfaces to allow model substitutions where desired.


### Simulation Flow <a name="simulation-flow"></a>

The GSM simulates particle-particle interactions by first simulating an intra-nuclear cascade via the standard or modified Dubna cascade models. Interaction progeny then may coalesce together to form light ions. GSM assumes that the residual nuclei from the initial particle-particle interaction then begin to deexcite via preequilibrium and evaporation processes. Several light ions may be emitted from the residuals during preequilibrium deexcitation. GSM simulates residual evaporation of the equilibrated residual nuclei by considering particle evaporation and fission of the compound nucleus via the Generalized Evaporation Model, GEM2. A Fermi break-up deexcitation is utilized in place of preequilibrium and evaporation emissions for sufficiently small residual nuclei.


### Software Architecture <a name="software-architecture"></a>

The GSM utilizes primarily Fortran2003 and Fortran2008 standards, however legacy Fortran code still remains within the Modified DCM sub-model. GSM internally utilizes several sub-libraries that each implement a piece of the model physics described in the [simulation flow](#simulation-flow). Each sub-library utilizes an object-oriented software approach. Dependencies are abstracted from each implementation to help provide a flexible, robust, and parallelizable event generator. The architecture of the GSM thus allows for simultaneous particle event simulations.


## Build Requirements <a name="build-requirements"></a>

The GSM was migrated to an object-oriented framework using primarily GNU Fortran with GCC-6.3.0. GSM should built with most version of GFortran, however compilation by other compilers may be limited. Help in this area would be greatly appreciated.

The following table shows the compilers that can successfully compile GSM:

| Compiler         | Version       | Note                                                  |
| :--------------- | :------------ | :---------------------------------------------------- |
| GCC              | 4.6.1         | Regression differences seen with GCC4.7.2 and earlier |
| Intel            | Not Tested    | |
| PGI              | Not Tested    | |
| LLVM             | Not Tested    | |
| Cygwin           | Not Tested    | |


## Compilation and Installation <a name="compilation"></a>

A GSM executable can be created by navigating to the primary GSM top-level directory with any of the compilers supported (see the [build requirements](#build-requirements) section for more information). GSM uses a relatively basic CMake build system, **requiring** out-of-source builds. Executing CMake against the GSM will then generate a Makefile which can be used to build any one of the GSM's sub-models or the GSM itself.

> Note that GSM utilizes a CMake build system. This build system has not yet been tested on non-linux operating systems.

The GSM may be installed in a directory called `build` by executing the following:

```bash
mkdir ./build
cd ./build
cmake <options> ../
make
```

GSM module files are all placed within the `modules` folder of the build path. Compiled libraries are stored in the `lib` folder, and the standalone GSM executable in the `bin` folder. The GSM further requires several data files to perform simulations. These files are provided in the `my_GSM` folder with the GSM executable for standalone simulations. The location of these files will likely need to be provided for simulations performed by consumers of the GSM upon data instantiation.


### Testing <a name="testing"></a>
Upon creation of the GSM executable, users should perform all provided regression tests to verify the proper compilation f the GSM event generator. The folder `test`, housed in the top-level GSM directory, contains all available tests for users.  
A script is provided to perform all regression tests to compare the results produced by GSM when built by `GCC6.3.0`. A limited number of tests are currently provided (30 at this time) with small event limits. It is recognized that the provided regression suite is not particularly thorough, however it is assumed to be sufficient to catch the majority of unintended changes.
 
Users can test their installation of GSM using the provided regression script from the top-level directory of GSM by typing the following:
```bash
sh -e test/regression/bin/regression.sh
```

The script will attempt to create a clean installation of GSM by creating a `./build/` directory, compiling GSM, and performing all regression tests.

> NOTE: The regression script utilizes a Python module (called `Comparer.py`) to compare output files from the previous the current versions of GSM. This python script will create a file for *each* set of input files tested stating the differences in the output file and flagging them. This Python script was developed utilizing `Python3.6`. It is possible that regression tests will successfully complete with an earlier version of Python, however it is not guaranteed to work properly.

Developers of the GSM will want to create a regression test baseline upon initial installation. Developers should ensure that the GSM builds successfully, and upon a successful build should execute the regression scripts and save the output files accordingly. See the regression bash file for details. Regression tests can be updated easily with the following command, again being typed in the top level directory of GSM, to overwrite the existing regression set with the last result set.
```bash
sh -e test/regression/bin/updateOutputs.sh
```


## Usage <a name="usage"></a>

Upon being compiled (see the [Compilation](#compilation) and [tested](#testing) sections), users, clients, or both can create and utilize GSM for their spallation simulations.


### Standalone (User) Usage <a name="standalone-usage"></a>

The GSM, as a standalone code, is contained by the `standaloneGSM` [module](src/DriverGSM/standaloneGSM.f90). An exectuble, named `xgsm?`, where the `?` is the current major version of GSM, is created upon successful [compilation](#compilation) in a directory named `my_GSM`, named accordingly in the [CMakeLists.txt file](./CMakeLists.txt).

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

No client usage of the GSM is presently recorded. Potential clients that may benefit from implementation of the GSM include, but are not limited to:
+ MCNP
+ GEANT
+ Others?


## Documentation <a name="docs"></a>

Documentation of the GSM prior to its creation may be found by referencing various CEM and LAQGSM event generator reports. Note that some of the documentation contained by GSM refers to some of these reports and the references therein. Some documentation of GSM developments may additionally be found from brief research. Documentation of the source code for GSM may be easily generated via Doxygen as follows:

```bash
cd ./Docs
doxygen Doxyfile.in
```

Doing so generates, by default, documentation in both a `LaTeX` and `HTML` format. The generated `HTML` documentation may be viewed by pointing a browser to the `index.html` file within the `Docs/Docs_HTML` directory.


## Contribute <a name="contribute"></a>

Contributions to GSM can be made by contacting Chase Juneau (ChaseJuneau@isu.edu) or Dr. Leslie Kerby (LeslieKerby@isu.edu) at the Idaho State University, or by contacting the [XCP-3](https://www.lanl.gov/org/padwp/adx/computational-physics/monte-carlo/index.php) group leader of the Los Alamos National Laboratory with information regarding the contribution.  

>Tip: developers can obtain a general understanding of the code styling within GSM by looking through the [Coalescence model code](./src/Coalescence). The code represented there is sufficiently short and clear to provide beginning developers with style reference and the architecture utilized for the GSM code.  Generally speaking, models are divided into (a) a parameter module, (b) a data module containing static data for a given simulation, and (c) a model implementation.

Questions about GSM, regarding code or physics, can be directed to C. Juneau (ChaseJuneau@isu.edu) and Dr. Leslie Kerby (LeslieKerby@isu.edu).


## Credits <a name="credits"></a>

GSM was derived from CEM03.03 and LAQGSM03.03. The primary authors of these event generators include Drs. S. Mashnik, K. Gudima, and A. Sierk, with important contributions from R. Prael, M. Baznat, and N. Mokhov.
Dr. L. Kerby, under the guidance of S. Mashnik, provided several important upgrades made to the CEM03.03 for light fragment emission in both the coalescense and preequilibrium implementations.
C. Juneau then combined the CEM03.03 and LAQGSM03.03 models to form the GSM under the guidance of L. Kerby. The GSM was then migrated to the current architecture by C. Juneau under the auspices of the Idaho State University and, in part, from funding provided by the Los Alamos National Laboritory under Idaho State Univeristy subcontract number 385443.

