# **The Fission Barrier Library**
---
The Fission Barrier module may be used by clients to obtain various ground state and fission barrier energies [MeV]. The Fission Barrier module can be used for fission barriers [MeV], rotating ground state energies [MeV], and fission barrier approximations given the atomic and mass numbers of a nucleus alongside the angular momentum, in units of h-bar.
The Fission Barrier library is part of the GSM project.  
`Code`: Fission Barrier utilizes strictly Fortran2003 and Fortran2008.  
`Parallel`: Fission Barrier is fully parallelizable due to its object-oriented implementation.

> **Fission Barrier Version 1.0.0**




## **Table of Contents**
___
1. [Disclaimer](#disclaimer)
2. [License](#license)
3. [Build Requirements](#build-requirements)
4. [Compilation](#compilation)
   -  [Testing](#testing)
5. [Usage](#usage)
   - [Client Usage](#client-usage)
6. [Contribute](#contribute)
7. [Credits](#credits)



## **Disclaimer <a name="disclaimer"></a>**
___
>NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL SECURITY LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.



## **License <a name="license"></a>**
___
The predecessor of GSM, CEM, has been approved for release with associated LA-CC number LA-CC-04-085. The GSM event generator has **not** been licensed at this time for public release.
The CEM material was procuded under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory for the U.S. Department of Energy. The Government is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this material to reproduce, prepare derivative works, and perform publicly and display publicly.  

**The GSM event generator code has not, at this time, been approved for public release under the CEM license. The Fission Barrier module is based on the GSM event generator.**  


## **Build Requirements <a name="build-requirements"></a>**
___
Fission Barrier was migrated to an object-oriented framework using primarily GNU Fortran using GCC-6.3.0. The following table shows the compilers that can successfully compile the GSM Fission Barrier library:

| Compiler         | Version       |
| :--------------- | :------------ |
| GCC              | 4.6.1         |
| Intel            | Not Tested    |
| PGI              | Not Tested    |
| LLVM             | Not Tested    |
| Cygwin           | Not Tested    |



## **Compilation <a name="compilation"></a>**
A Fission Barrier library can be created by navigating to the primary GSM top-level directory with any of the compilers supported (see the [build requirements](#build-requirements) section for more information). GSM uses a relatively basic CMake
 build system, **requiring** out-of-source builds. To compile the Fission Barrier sub-library, type the following in a terminal:
```
mkdir <build-src>
cd <build-src>
cmake <options> <path-to-GSM-top-level-directory>
make gsm_fissionBarrier
```

More explicitly stated, users type the following to build the Fission Barrier library in the top-level directory of GSM:
```
mkdir ./build
cd ./build
cmake ..
make gsm_fissionBarrier
```

> Note: the CMake build platform was created on a Linux operating system; this CMake distribution has NOT been tested on non-linux operating systems.


Upon compilation of the Fission Barrier library, all associated `*.mod` files are stored in the build directory's `./modules` folder, and the Fission Barrier library object is stored in the `build` directory's `./lib` folder. Note that the following modules are part of the Fission Barrier library:   
- fissionBarrierParams.f90
- fissionBarrierClass.f90 

> Note: The fission barrier object `USE`s the GSM Molnix library. Building the Fission Barrier library will also build the Molnix library.

### **Testing <a name="testing"></a>**
No testing is provided for by the Fission Barrier object or its modules. Proper compilation currently is only tested with the totality of GSM using the provided regression suite. Note that this regression suite does *not* provide any details regarding potential differences in numerics regarding the Fission Barrier object.  
Upon creation of the GSM executable (see the [Compilation](#compilation) section), users should perform all provided regression tests to verify the proper compilation of the GSM event generator.



## **Usage <a name="usage"></a>**
---
Upon being compiled and tested (see the [Compilation](#compilation) and [testing](#testing) sections), clients may create and utilize the Fission Barrier class for their simulations by USEing the `fissionBarrierClass` module and linking, during compilation, against the `gsm_molnix` library.


### **Client Usage (Fission Barrier API) <a name="client-usage"></a>**
An instance of the Fission Barrier class is created by USEing the `FissionBarrier` object and the `newFissionBarrier` constructor from the `fissionBarrierClass` module. Construction of the fission barrier object is required. Users must pass a `Molnix` object in to the `newFissionBarrier` constructor that the fission barrier object will point to. Clients may optionally pass in a `fissionBarrierOptions` data type, from the same module, to control the behavior of the initialized Fission Barrier object. Clients to the Fission Barrier object may also control all messages produced from the Fission Barrier object by passing in a procedure pointer with an appropriate interface (see the sample API below). Options passed in through this object are validated upon object construction. The sample API below provides an example of construction and usage of the fission barrier object:  

```fortran
subroutine sampleFissionBarrierAPI( molnixObj )
   ! Import fission barrier class and appropriate members:
   use molnixClass, only: Molnix
   use fissionBarrierClass, only: &
         & FissionBarrier,        &   ! The object
         & newFissionBarrier,     &   ! The constructor
         & fissionBarrierOptions      ! Fission Barrier options

   implicit none
   type(Molnix), intent(in   ) :: molnixObj

   ! Create interface to create procedure pointer with:
   abstract interface
      subroutine IOHANDLER(msgVerbose, msgType, msg)
         use, intrinsic:: iso_fortran_env, only: int32
         integer(int32),   intent(in) :: msgVerbose   ! Sets importance of message
         integer(int32),   intent(in) :: msgType      ! Type of message (fatal error, error, warning, or comment)
         character(len=*), intent(in) :: msg          ! A string containing the message
      end subroutine IOHANDLER
   end interface
   procedure(IOHANDLER), pointer :: newPrint => clientMsgHandler   ! OPTIONAL constructor argument
   type(fissionBarrierOptions) :: options   ! OPTIONAL constructor argument


   ! Construct Fission Barrier object with all optional arguments
   ! type(FissionBarrier) :: theFissBarObject = newFissionBarrier( molnixObj )   ! Minimum required arguments for construction
   type(FissionBarrier) :: theFissBarObject = newFissionBarrier( molnixObj, options, newPrint )

   return
end subroutine sampleFissionBarrierAPI
```

The Fission Barrier class has several accessible procedures that clients may utilize. The following procedures that may be used to obtain information regarding a nucleus's rotating ground state energies [MeV] and fission barrier energies [MeV] include:

| Procedure    | Arguments                                      | Description                                                                          |
| :----------- | :--------------------------------------------- | :----------------------------------------------------------------------------------- |
| bf           | A, Z [real64], ln [int32], egs0 [real64]       | Returns the fission barrier energy, in [MeV], in the `egs0` argument                 |
| barfit       | A, Z [real64], ln [int32], bfis, egs [real64]  | Returns the fission barrier energy and rotating ground state energies [MeV] in the `bfis` and `egs` arguments, respectively |
| bsfit        | Z, A [real64], ln [int32], bs [real64]         | Returns the saddle-point fission barrier energy [MeV], in the `bs` argument, from an evaluated 3-D Bs(Z,A,L) fitting function |

The remainder of the procedures clients may use to query the Fission Barrier object and its data include:  

| Procedure           | Arguments      | Description                                                                      |
| :------------------ | :------------- | :------------------------------------------------------------------------------- |
| properlyConstructed | None           | Returns logical flag stating if the object was constructed                       |
| queryOptions        | None           | Returns the `fissionBarrierOptions` data type held by the Fission Barrier object |
| queryMolnix         | None           | Returns the `Molnix` object pointer held by the Fission Barrier object           |


The variable members contained within the `fissionBarrierOptions` data type include the following:

| Variable Name   | Type     | Default   | Description                                                |
| :-------------- | :------- | :-------- | :--------------------------------------------------------- |
| r0m             | real64   | 1.20      | Scaling multiplier for determining the nucleus's radius (r = r0 * A^{1/3}) |


## **Contribute <a name="contribute"></a>**
---
Contributions to GSM can be made by contacting Dr. Leslie Kerby (kerblesl@isu.edu) at the Idaho State University of by contacting the [XCP-3](https://www.lanl.gov/org/padwp/adx/computational-physics/monte-carlo/index.php) group leader, currently Avneet Sood (sooda@lanl.gov), of the Los Alamos National Laboratory with information regarding the contribution.  

Questions about GSM, regarding code or physics, can be directed to C. Juneau (junechas@isu.edu) and Dr. Leslie Kerby (kerblesl@isu.edu).


## **Credits <a name="credits"></a>**
---
The primary authors of GSM, originally CEM, include Drs. S. Mashnik, K. Gudima, and A. Sierk, with important contributions from R. Prael, M. Baznat, and N. Mokhov. Upgrades made to CEM03.03 prior to GSM's creation can be mostly accredited to Dr. Leslie Kerby, who provided several important updates to the preequilibrium model of GSM. GSM was migrated to an object-oriented framework by C. Juneau under the auspices of the Idaho State University and, in part, from funding provided by the Los Alamos National Laboritory under Idaho State University subcontract number 385443.