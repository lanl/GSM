# **The Coalescence Library**
---
The Coalescence module may be used by clients to simulate the compounding of secondary cascade nucleons. The Coalescence object utilizes coalescence radii to approximate the quantum wave overlap for various secondary cascade nucleons for their potential to coalesce.
The Coalescence library is part of the GSM project.  
`Code`: Coalescence utilizes strictly Fortran2003 and Fortran2008.  
`Parallel`: Coalescence is fully parallelizable due to its object-oriented implementation.

> **Coalescence Version 1.0.0**




## **Table of Contents**
___
1. [Disclaimer](#disclaimer)
2. [License](#license)
3. [Build Requirements](#build-requirements)
4. [Compilation](#compilation)
   -  [Testing](#testing)
5. [Usage](#usage)
   - [Client Usage](#client-usage)
   - [Data Types](#dataTypes)
6. [Contribute](#contribute)
7. [Credits](#credits)



## **Disclaimer <a name="disclaimer"></a>**
___
>NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL SECURITY LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.



## **License <a name="license"></a>**
___
The predecessor of GSM, CEM, has been approved for release with associated LA-CC number LA-CC-04-085. The GSM event generator has **not** been licensed at this time for public release.
The CEM material was procuded under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory for the U.S. Department of Energy. The Government is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this material to reproduce, prepare derivative works, and perform publicly and display publicly.  

**The GSM event generator code has not, at this time, been approved for public release under the CEM license. The Coalescence module is based on the GSM event generator.**  


## **Build Requirements <a name="build-requirements"></a>**
___
Coalescence was migrated to an object-oriented framework using primarily GNU Fortran using GCC-6.3.0. The following table shows the compilers that can successfully compile the GSM Coalescence library:

| Compiler         | Version       |
| :--------------- | :------------ |
| GCC              | 4.6.1         |
| Intel            | Not Tested    |
| PGI              | Not Tested    |
| LLVM             | Not Tested    |
| Cygwin           | Not Tested    |



## **Compilation <a name="compilation"></a>**
A Coalescence library can be created by navigating to the primary GSM top-level directory with any of the compilers supported (see the [build requirements](#build-requirements) section for more information). GSM uses a relatively basic CMake build system, **requiring** out-of-source builds. To compile the Coalescence sub-library, type the following in a terminal:
```
mkdir <build-src>
cd <build-src>
cmake <options> <path-to-GSM-top-level-directory>
make gsm_coalescence
```

More explicitly stated, users type the following to build the Coalescence library in the top-level directory of GSM:
```
mkdir ./build
cd ./build
cmake ..
make gsm_coalescence
```

> Note: the CMake build platform was created on a Linux operating system; this CMake distribution has NOT been tested on non-linux operating systems.


Upon compilation of the Coalescence library, all associated `*.mod` files are stored in the build directory's `./modules` folder, and the Coalescence library object is stored in the `build` directory's `./lib` folder. Note that the following modules are part of the Coalescence library:   
- coalescenceParams.f90
- coalescenceClass.f90 


### **Testing <a name="testing"></a>**
No testing is provided for by the Coalescence object or its modules. Proper compilation currently is only tested with the totality of GSM using the provided regression suite. Note that this regression suite does *not* provide any details regarding potential differences in numerics regarding the Coalescence object.  
Upon creation of the GSM executable (see the [Compilation](#compilation) section), users should perform all provided regression tests to verify the proper compilation of the GSM event generator.



## **Usage <a name="usage"></a>**
---
Upon being compiled and tested (see the [Compilation](#compilation) and [testing](#testing) sections), clients may create and utilize the Coalescence class for their simulations by USEing the `coalescenceClass` module and linking to it during compilation.


### **Client Usage (Coalescence API) <a name="client-usage"></a>**
An instance of the Coalescence class is created by USEing the `Coalescence` object and the `newCoalescence` constructor from the `coalescenceClass` module. Construction of the coalescence object is required. Users are highly encouraged to pass a `CoalescenceData` object in to the `newCoalescence` constructor that the coalescence object will copy for its simulations.

> The `CoalescenceData` object may be constructed by passing in the incident energy of the reaction, in [GeV], to establish default coalescence radii. Clients may optionally pass in the `coalescenceDataOptions` data type to specify their own coalescence radii, if [>0]. The `CoalescenceData` object utilizes default coalescence radii when not specified by the client.

Clients may also optionally pass in a `coalescenceOptions` data type, from the same module, to control the behavior of the initialized Coalescence object. Clients to the Coalescence object may also control all messages produced from the Coalescence object by passing in a procedure pointer with an appropriate interface (see the sample API below). Options passed in through this object are validated upon object construction. The sample API below provides an example of construction and usage of the coalescence data object and simulation object:  

```fortran
subroutine sampleCoalescenceDATAAPI( incidentEnergyGeV )

   ! Import coalescence data objects:
   use, intrinsic:: iso_fortran_env, only: real64
   use coalescenceClass, only: &
        & CoalescenceData,    &    ! Data object
	& newCoalescenceData, &    ! Constructor for data object
	& coalescenceDataOptions   ! Data object options

   implicit none
   real(real64), intent(in   ) :: incidentEnergyGeV

   ! This object may be used to control the various coalescence radii used -
   ! NOTE: Negative values are ignored for this
   type(coalescenceDataOptions) :: options   ! OPTIONAL constructor argument


   ! Create objects:
   ! type(CoalescenceData) :: coalesDataObj = newCoalescenceData( incidentEnergyGeV )   ! Minimum required arguments for construction
   type(CoalescenceData) :: coalesDataObj = newCoalescenceData( incidentEnergyGeV, options )


   return
end subroutine sampleCoalescenceDATAAPI
```

```fortran
subroutine sampleCoalescenceAPI( coalesDataObj )

   ! Import coalescence class and appropriate members:
   use, intrinsic:: iso_fortran_env, only: real64
   use coalescenceClass, only: &
         & CoalescenceData,    &   ! The data object
         & Coalescence,        &   ! The object
         & newCoalescence,     &   ! The constructor
         & coalescenceOptions      ! Coalescence options

   implicit none
   type(CoalescenceData), intent(in   ) :: coalesDataObj   ! OPTIONAL constructor argument

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
   type(coalescenceOptions) :: options   ! OPTIONAL constructor argument


   ! Construct Coalescence object with all optional arguments
   ! type(Coalescence) :: theCoalesObject = newCoalescence( )   ! Minimum required arguments for construction
   type(Coalescence) :: theCoalesObject = newCoalescence( coalesDataObj, options, newPrint )

   return
end subroutine sampleCoalescenceyAPI
```

The Coalescence data class has several accessible procedures that clients, such as the Coalescence class, may utilize. These procedures are used primarily to query the coalescence radii stored within the data object, however others exist as well. These procedures include the following:

| Procedure               | Arguments     | Description                                                       |
| :---------------------- | :------------ | :---------------------------------------------------------------- |
| coalesradiiDeut         | None          | Returns the coalescence radius for compounding deuterons          |
| coalesradiiTrit         | None          | Returns the coalescence radius for compounding tritons and helion |
| coalesradiiAlpha        | None          | Returns the coalescence radius for compounding alphas             |
| coalesradiiLFrag        | None          | Returns the coalescence radius for compounding light fragments    |
| properlyConstructedData | None          | Returns a logical flag stating construction state of the object   |
| properlyConstructed     | Same as `properlyConstructedData` | |

The Coalescence class has several accessible procedures that clients may utilize. The following procedures that may be used include:

| Procedure           | Arguments                    | Description                                                   |
| :------------------ | :--------------------------- | :------------------------------------------------------------ |
| simulateCoalescence | results [coalescenceResults] | Simulates coalescence of progeny in the `results` data type   |
| simulate            | Same as `simulateCoalescence`| |
| execute             | Same as `simulateCoalescence`| |
| start               | Same as `simulateCoalescence`| |
| coalesce            | Same as `simulateCoalescence`| |
| compound            | Same as `simulateCoalescence`| |

A coalescence simulation is completed by first declaring and constructing a `coalescenceResults` object. All progeny to be considered for coalescence should be loaded into the `coalescenceResults` object prior to simulation as the coalescence object looks at loaded progeny for its simulation. The sample usage API below demonstrates the declaration, construction, and usage of the `coalescenceResults` object.

```fortran
subroutine sampleCoalescenceRESULTSAPI( coalesObj )

   use coalescenceClass, only: &
        & Coalescence,         &   ! Primary object
	& coalescenceParticle, &   ! Particle object
	& coalescenceResults,  &   ! Results object for Coalescence simulation
	& newCoalescenceResults 

   implicit none
   type(Coalescence), intent(inout) :: coalesObj


   ! Create an array of particles:
   type(coalescenceParticle), dimension( numINCProgeny ) :: particleBnk

   ! Setup results object:
   type(coalescenceResults) :: results = newCoalescenceResults( particleBnk )

   ! -----------------------------------------------------------------------------------------
   ! NOTE: At this point in the API, clients would want to store particle information into the
   !       `particleBnk` or `results%partBnk` arrays [the results object points to the passed
   !       in array] from their own arrays of particle information
   ! -----------------------------------------------------------------------------------------

   ! Perform coalescence simulation:
   call coalesObj%coalesce( results )

   ! -----------------------------------------------------------------------------------------
   ! NOTE: At this point in the API, clients would want to store particle information from the
   !       `particleBnk` or `results%partBnk` arrays [the results object points to the passed
   !       in array] into their own particle tracking arrays, if present.
   ! -----------------------------------------------------------------------------------------

   return
end subroutine sampleCoalesenceRESULTSAPI
```

The remainder of the procedures clients may use to query the Coalescence object and its data include:  

| Procedure                 | Arguments      | Description                                                           |
| :------------------------ | :------------- | :-------------------------------------------------------------------- |
| properlyConstructedObject | None           | Returns logical flag stating if the object was constructed            |
| properlyConstructed       | Same as `properlyConstructedObject` | |
| queryOptions              | None           | Returns the `coalescenceOptions` data type held by the Coalescence object   |
| queryData                 | none           | Returns the `CoalescenceData` class held by the Coalescence object    |


### **Data Types <a name="dataTypes"></a>**
The coalescence class utilizes several data types for its simulations. These data types include `coalescenceOptions`, `coalescenceProgeny`, and `coalescenceResults`.   

The variable members contained within the `coalescenceDataOptions` data type include the following:

| Variable Name     | Type         | Default   | Description                                                   |
| :---------------- | :----------- | :-------- | :------------------------------------------------------------ |
| coalesRadiiDeut   | real(real64) | 0.090     | Coalescence radius for compounding deuterons [GeV/c]          |
| coalesRadiiTrit   | real(real64) | 0.108     | Coalescence radius for compounding tritium and helion [GeV/c] |
| coalesRadiiAlpha  | real(real64) | 0.130     | Coalescence radius for compounding alphas [GeV/c]             |
| coalesRadiiLFrag  | real(real64) | 0.175     | Coalescence radius for compounding light fragments [GeV/c]    |


The variable members contained within the `coalescenceOptions` data type include the following:

| Variable Name       | Type     | Default   | Description                                                                  |
| :------------------ | :------- | :-------- | :--------------------------------------------------------------------------- |
| expandedCoalescence | integer(int32)    | 2         | Flags use (>0) of the expanded coalescence model to compound light fragments |


The variable members contained within the `coalescenceParticle` data type include the following:  

| Variable Name          | Type           | Description                                       |
| :--------------------- | :------------- | :------------------------------------------------ |
| numBaryons             | integer(int32) | Number of baryons in the particle                 |
| charge                 | integer(int32) | Number of protons (i.e. charge) in the particle   |
| strangeness            | integer(int32) | Quantum decay number of the particle              |
| kinEnergy              | rea(real64)    | Kinetic energy [GeV] of the particle              |
| restMass               | rea(real64)    | Rest mass [GeV/c^2] of the particle               |
| linearMomX             | rea(real64)    | X-component of the particle's linear momentum     |
| linearMomY             | rea(real64)    | Y-component of the particle's linear momentum     |
| linearMomZ             | rea(real64)    | Z-component of the particle's linear momentum     |
| coalesceFlag           | rea(real64)    | Flags a fragment that was formed from coalescence |
| coalesceNum  [private] | integer(int32) | Flags pairs of particles that are to coalesce     |


The variable members contained within the `coalescenceParticle` data type include the following:  

| Variable Name          | Type                | Description                                                 |
| :--------------------- | :------------------ | :---------------------------------------------------------- |
| constructed  [private] | logical             | Flags if the object was contructed by the client            |
| simState               | integer(int32)      | Flags how the simulation ended (=0 for no warnings)         |
| partBnk                | coalescenceParticle | Array of particles used for the coalescence simulation      |
| numParticles           | integer(int32)      | Number of particles in the `partBnk` array                  |
| partBnkSize            | integer(int32)      | Maximum number of particles that may exist                  |
| numCoalesced           | integer(int32)      | States the number of fragments that coalesced in simulation |
| numFormedFragments     | integer(int32)      | States the number of compounded fragments in simulation     |



## **Contribute <a name="contribute"></a>**
---
Contributions to GSM can be made by contacting Dr. Leslie Kerby (kerblesl@isu.edu) at the Idaho State University of by contacting the [XCP-3](https://www.lanl.gov/org/padwp/adx/computational-physics/monte-carlo/index.php) group leader, currently Avneet Sood (sooda@lanl.gov), of the Los Alamos National Laboratory with information regarding the contribution.  

Questions about GSM, regarding code or physics, can be directed to C. Juneau (junechas@isu.edu) and Dr. Leslie Kerby (kerblesl@isu.edu).


## **Credits <a name="credits"></a>**
---
The primary authors of GSM, originally CEM, include Drs. S. Mashnik, K. Gudima, and A. Sierk, with important contributions from R. Prael, M. Baznat, and N. Mokhov. Upgrades made to CEM03.03 prior to GSM's creation can be mostly accredited to Dr. Leslie Kerby, who provided several important updates to the preequilibrium model of GSM. GSM was migrated to an object-oriented framework by C. Juneau under the auspices of the Idaho State University and, in part, from funding provided by the Los Alamos National Laboritory under Idaho State University subcontract number 385443.