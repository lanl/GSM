# **The Molnix Library**
---
The Molnix module may be used by clients to obtain various energy values. The supported energy values include mass excesses, shell energies, and pairing gap energies. Molnix is part of the GSM project.  
`Code`: Molnix utilizes strictly Fortran2003 and Fortran2008.  
`Parallel`: Molnix is fully parallelizable due to its object-oriented implementation.

> **Molnix Version 1.0.0**




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

**The GSM event generator code has not, at this time, been approved for public release under the CEM license. The Molnix module is based on the GSM event generator.**  


## **Build Requirements <a name="build-requirements"></a>**
___
Molnix was migrated to an object-oriented framework using primarily GNU Fortran using GCC-6.3.0. The following table shows the compilers that can successfully compile the GSM Molnix library:

| Compiler         | Version       |
| :--------------- | :------------ |
| GCC              | 4.6.1         |
| Intel            | Not Tested    |
| PGI              | Not Tested    |
| LLVM             | Not Tested    |
| Cygwin           | Not Tested    |



## **Compilation <a name="compilation"></a>**
A Molnix library can be created by navigating to the primary GSM top-level directory with any of the compilers supported (see the [build requirements](#build-requirements) section for more information). GSM uses a relatively basic CMake
 build system, **requiring** out-of-source builds. To compile the Molnix sub-library, type the following in a terminal:
```
mkdir <build-src>
cd <build-src>
cmake <options> <path-to-GSM-top-level-directory>
make gsm_molnix
```

More explicitly stated, users type the following to build the Molnix library in the top-level directory of GSM:
```
mkdir ./build
cd ./build
cmake ..
make gsm_molnix
```

> Note: the CMake build platform was created using Linux; this CMake distribution has NOT been tested on non-linux operating systems.


Upon compilation of the Molnix library, all associated `*.mod` files are stored in the build directory's `./modules` folder, and the Molnix library object is stored in the `build` directory's `./lib` folder. Note that the following modules are part of the Molnix library:   
- molnixParams.f90
- molnixClass.f90 


### **Testing <a name="testing"></a>**
No testing is provided for by the Molnix object or its modules. Proper compilation currently is only tested with the totality of GSM using the provided regression suite. Note that this regression suite does *not* provide any details regarding potential differences in numerics regarding the Molnix object.  
Upon creation of the GSM executable (see the [Compilation](#compilation) section), users should perform all provided regression tests to verify the proper compilation of the GSM event generator.



## **Usage <a name="usage"></a>**
---
Upon being compiled and tested (see the [Compilation](#compilation) and [testing](#testing) sections), clients may create and utilize the Molnix class for their simulations by USEing the `molnixClass` module and linking, during compilation, against the `gsm_molnix` library.


### **Client Usage (Molnix API) <a name="client-usage"></a>**
An instance of the Molnix class is created by USEing the `Molnix` object and the `newMolnix` constructor from the `molnixClass` module. Clients may optionally pass in a `molnixOptions` data type, from the same module, to control the behavior of the initialized Molnix object. Options passed in through this object are validated upon object construction. The Molnix object does not *require* construction, however it is recommended as clients may take control of the molnix object's message printing procedure pointer.  

```fortran
subroutine sampleMolnixAPI()

   ! Import molnix class and appropriate members:
   use molnixClass, only: &
         & Molnix,        &   ! The object
         & newMolnix,     &   ! The constructor
         & molnixOptions      ! Molnix options

   implicit none

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
   type(molnixOptions) :: options   ! OPTIONAL constructor argument


   ! Construct Molnix object with all optional arguments
   ! type(Molnix) :: theMolnixObject = newMolnix( )   ! Minimum required arguments for construction
   type(Molnix) :: theMolnixObject = newMolnix( options, newPrint )

   return
end subroutine sampleMolnixAPI
```

The Molnix class has several accessible procedures that clients may utilize. The following procedures that may be used to obtain information regarding a nucleus's  mass excess, shell energies, or pairing gaps include:

| Procedure    | Arguments           | Description                                                                          |
| :----------- | :------------------ | :----------------------------------------------------------------------------------- |
| massExcess   | A, Z [real64]       | Returns the experimental, or theoretical if not measure, mass excess [MeV]           |
| shellEnergy  | A, Z [real64]       | Returns ground state shell and pairing corrections [MeV]                             |
| pairingGap   | A, Z [real64]       | Returns the energy shift from the Moller, Nix, & Kratz calculated pairing gaps [MeV] |
| defineEnergy | Z, N, type [int32]  | Returns similarly to shellEnergy, massExcess, or pairingGap for `type=1`, `2`, or `3`, respectively [MeV] |

The remainder of the procedures clients may use to query the Molnix object and its data include:  

| Procedure           | Arguments      | Description                                                |
| :------------------ | :------------- | :--------------------------------------------------------- |
| properlyConstructed | None           | Returns logical flag stating if the object was constructed |
| queryOptions        | None           | Returns the `molnixOptions` type held by the Molnix object |
| nmina               | index [int32]  | Returns the value of `nmina` parameterized array           |
| nmaxa               | index [int32]  | Returns the value of `nmaxa` parameterized array           |
| nmin                | index [int32]  | Returns the value of `nmin`  parameterized array           |
| nmax                | index [int32]  | Returns the value of `nmax`  parameterized array           |


The variable members contained within the `molnixOptions` data type include the following:

| Variable Name   | Type     | Default   | Description                                                |
| :-------------- | :------- | :-------- | :--------------------------------------------------------- |
| cevap           | real64   | 12.0      | Pairing energy gap scaling factor for theoretical values   |


## **Contribute <a name="contribute"></a>**
---
Contributions to GSM can be made by contacting Dr. Leslie Kerby (kerblesl@isu.edu) at the Idaho State University of by contacting the [XCP-3](https://www.lanl.gov/org/padwp/adx/computational-physics/monte-carlo/index.php) group leader, currently Avneet Sood (sooda@lanl.gov), of the Los Alamos National Laboratory with information regarding the contribution.  

Questions about GSM, regarding code or physics, can be directed to C. Juneau (junechas@isu.edu) and Dr. Leslie Kerby (kerblesl@isu.edu).


## **Credits <a name="credits"></a>**
---
The primary authors of GSM, originally CEM, include Drs. S. Mashnik, K. Gudima, and A. Sierk, with important contributions from R. Prael, M. Baznat, and N. Mokhov. Upgrades made to CEM03.03 prior to GSM's creation can be mostly accredited to Dr. Leslie Kerby, who provided several important updates to the preequilibrium model of GSM. GSM was migrated to an object-oriented framework by C. Juneau under the auspices of the Idaho State University and, in part, from funding provided by the Los Alamos National Laboritory under Idaho State University subcontract number 385443.