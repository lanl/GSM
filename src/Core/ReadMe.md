# Core

![Contributors](https://img.shields.io/github/contributors/lanl/generalized-spallation-model)
![Stars](https://img.shields.io/github/stars/lanl/generalized-spallation-model)
![License](https://img.shields.io/badge/license-BSD3-green)
![Forks](https://img.shields.io/github/forks/lanl/generalized-spallation-model)
![Issues](https://img.shields.io/github/issues/lanl/generalized-spallation-model)


The core library provides a functional set of simple and robust tools upon which other libraries may base themselves. The abstracted approach provides for flexibility within consuming libraries.

> Core is not meant to be completely comprehensive, but acts as a foundation upon which other libraries may be built. It should contain ground-level tools that can be used throughout a consuming library.


## Core Consumption

The core library may be consumed by other libraries by `use`-ing the module, e.g.:

```Fortran
use GSM_Core

implicit none
real(gsmDouble), parameter :: example = zro
...
```


### Core Modules

| Module    | Description                             |
| :-------- | :-------------------------------------- |
| Types     | Defines generic Fortran and C-style numeric values |
| Numbers   | Defines parameterized numbers, physical constants, and other generic parameterizations |
| Macro     | Light-weight macro functionality (e.g., limited math protection, debugging, etc.) |
| DBC       | Design-by-Contract module to ensure developmental requirements are satisfied during the development process. Not typically applicable in production. |
| FException | Abstracted exception handling wrapper for clients to optionally implement if desired. |


----------------------------------------------------------

## Contribution Guidelines

Please review the project's [contrbiution guidelines](../../Contribue.md] if you'd like to contribute to the project.
