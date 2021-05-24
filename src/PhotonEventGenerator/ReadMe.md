# **The Photon Event Generator Library**
---
The Photon Event Generator module may be used by clients to determine cross section information for various photon reactions.
The Photon Event Generator module can be used to interpolate between angular photon cross section distributions, obtain photon cross sections for gamma+A reactions, and perform Lorentz transformations.
The Photon Event Generator library is part of the GSM project.  
`Code`: Photon Event Generator utilizes strictly Fortran2003 and Fortran2008.  
`Parallel`: Photon Event Generator is fully parallelizable due to its object-oriented implementation.

> **Photon Event Generator Version 1.0.0**




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

**The GSM event generator code has not, at this time, been approved for public release under the CEM license. The Photon Event Generator module is based on the GSM event generator.**  


## **Build Requirements <a name="build-requirements"></a>**
___
Photon Event Generator was migrated to an object-oriented framework using primarily GNU Fortran using GCC-6.3.0. The following table shows the compilers that can successfully compile the GSM Photon Event Generator library:

| Compiler         | Version       |
| :--------------- | :------------ |
| GCC              | 4.6.1         |
| Intel            | Not Tested    |
| PGI              | Not Tested    |
| LLVM             | Not Tested    |
| Cygwin           | Not Tested    |



## **Compilation <a name="compilation"></a>**
A Photon Event Generator library can be created by navigating to the primary GSM top-level directory with any of the compilers supported (see the [build requirements](#build-requirements) section for more information). GSM uses a relatively basic CMake
 build system, **requiring** out-of-source builds. To compile the Photon Event Generator sub-library, type the following in a terminal:
```
mkdir <build-src>
cd <build-src>
cmake <options> <path-to-GSM-top-level-directory>
make gsm_photonEventGenerator
```

More explicitly stated, users type the following to build the Photon Event Generator library in the top-level directory of GSM:
```
mkdir ./build
cd ./build
cmake ..
make gsm_photonEventGenerator
```

> Note: the CMake build platform was created on a Linux operating system; this CMake distribution has NOT been tested on non-linux operating systems.


Upon compilation of the Photon Event Generator library, all associated `*.mod` files are stored in the build directory's `./modules` folder, and the Photon Event Generator library object is stored in the `build` directory's `./lib` folder. Note that the following modules are part of the Photon Event Generator library:   
- photonEventGenerator.f90 


### **Testing <a name="testing"></a>**
No testing is provided for by the Photon Event Generator object or its modules. Proper compilation currently is only tested with the totality of GSM using the provided regression suite. Note that this regression suite does *not* provide any details regarding potential differences in numerics regarding the Photon Event Generator object.  
Upon creation of the GSM executable (see the [Compilation](#compilation) section), users should perform all provided regression tests to verify the proper compilation of the GSM event generator.



## **Usage <a name="usage"></a>**
---
Upon being compiled and tested (the see [Compilation](#compilation) and [testing](#testing) sections), clients may create and utilize the Photon Event Generator class for their simulations by USEing the `photonEventGeneratorClass` module.


### **Client Usage (Photon Event Generator API) <a name="client-usage"></a>**
An instance of the Photon Event Generator class is created by USEing the `PhotonEventGenerator` object and the `newPhotonEventGenerator` constructor from the `photonEventGeneratorClass` module. Construction of the photon event generator object is required. Users must pass an array of theta values, as well as an the integer size of the array, in to the `newPhotonEventGenerator` constructor that the photon event generator will use for interpolating angular distributions of photon cross sections.
Clients to the Photon Event Generator object may control all messages produced from the Photon Event Generator object by passing in a procedure pointer with an appropriate interface (see the sample API below). The sample API below provides an example of construction and usage of the photon event generator object:  

```fortran
subroutine samplePhotonEventGeneratorAPI( thetaValues, numElements )

   ! Import photon event generator class and appropriate members:
   use, intrinsic:: iso_fortran_env, only: int32, real64
   use photonEventGeneratorClass, only: &
         & PhotonEventGenerator,        &   ! The object
         & newPhotonEventGenerator          ! The constructor

   implicit none
   real(real64),   intent(in   ) :: thetaValues(:)
   integer(int32), intent(in   ) :: numElements


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


   ! Construct Photon Event Generator object with all optional arguments
   ! type(PhotonEventGenerator) :: thePhotoEGObject = newPhotonEventGenerator( thetaValues, numElements )   ! Minimum required arguments for construction
   type(PhotonEventGenerator) :: thePhotoEGObject = newPhotonEventGenerator( thetaValues, numElements, newPrint )

   return
end subroutine samplePhotonEventGeneratorAPI
```

The Photon Event Generator class has several accessible procedures that clients may utilize. The following procedures that may be used include:

| Procedure         | Arguments                              | Description                                                                                                |
| :---------------- | :------------------------------------- | :--------------------------------------------------------------------------------------------------------- |
| lorentzTransform  | p(3), v(3), pstar(3), t, cm [real64]   | Momentum calculation, stored in the `pstar` argument, in the system which has relative velocity `v` to given one |
| photoCrossSection | e, a0 [real64]                         | Returns the gamma+A cross section                                                                          |
| qintxs            | x, sig(m,n) [real64], l, m, n [int32]  | Returns the interpolated value of gamma+N differential cross sections                                      |


The remainder of the procedures clients may use to query the Photon Event Generator object and its data include:  

| Procedure           | Arguments      | Description                                                                      |
| :------------------ | :------------- | :------------------------------------------------------------------------------- |
| properlyConstructed | None           | Returns logical flag stating if the object was constructed                       |




## **Contribute <a name="contribute"></a>**
---
Contributions to GSM can be made by contacting Dr. Leslie Kerby (kerblesl@isu.edu) at the Idaho State University of by contacting the [XCP-3](https://www.lanl.gov/org/padwp/adx/computational-physics/monte-carlo/index.php) group leader, currently Avneet Sood (sooda@lanl.gov), of the Los Alamos National Laboratory with information regarding the contribution.  

Questions about GSM, regarding code or physics, can be directed to C. Juneau (junechas@isu.edu) and Dr. Leslie Kerby (kerblesl@isu.edu).


## **Credits <a name="credits"></a>**
---
The primary authors of GSM, originally CEM, include Drs. S. Mashnik, K. Gudima, and A. Sierk, with important contributions from R. Prael, M. Baznat, and N. Mokhov. Upgrades made to CEM03.03 prior to GSM's creation can be mostly accredited to Dr. Leslie Kerby, who provided several important updates to the preequilibrium model of GSM. GSM was migrated to an object-oriented framework by C. Juneau under the auspices of the Idaho State University and, in part, from funding provided by the Los Alamos National Laboritory under Idaho State University subcontract number 385443.