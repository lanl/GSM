# ExceptionLite

Provides a basic method for throwing errors to centralize error creation, management, and checking. This allows consumers to check for errors, proactively, without halting program execution in the event of an error.


## Description

The `ExceptionLite` module provided here creates a basic interface to create and manage errors. Doing so allows the following:
+ Centrally create errors (e.g., `call throw_insist_error("CONDITION", "MESSAGE", __FILE__, __LINE__)`)
+ Centrally look up errors and adjust behaviors at their presence (e.g., `if errors_thrown() then error stop "Errors were thrown" end if`)
+ Easily create new error types and create errors at runtime
+ Easily clear errors (should be done sparingly)
+ Provides simple, but specialized, errors for features that are not configured or implemented.

All lite exceptions should send the error messages to `error_unit` to inform the end-user of the error and why failure occurred. Additional details should be provided if the software was compiled under a non-release build.

