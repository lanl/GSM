# Contracts

Provides a basic method for enforcing conditions based on compilation level.


## Description

The `Contracts` provided here act as interface contracts to help ensure the validity of conditions. Such conditions can include:
+ Input validation
+ Output validation
+ Test Enforcement
+ Consistent stop behavior per Contract
+ Compilation-based enforcement (e.g., turn off some `contracts` based on compilation options)

The provided `contracts` additionally allow sending messages to the `error_unit` to inform the end-user of the error and why failure occurred. This is often useful for testing, but has many other use cases too.
