# GSM Regression Testing

This holds a base regression suite and bash scripts to aid in regression testing. Regression testing is most beneficial to developers to ensure any software changes do not break the GSM nor do they unexpectedly modify the numerical precision of simulated events.


## Purpose

The regression suite provided here is intended to help simplify deployment and development of the GSM. It provides a quick way for users and consumers to verify that any changes made to GSM do not break the GSM or modify expected results.

This regression suite utilizes bash and python scripts to help in simplifying simulation execution, comparing result sets, and reporting to the user.


## Execution

The regression suite may be executed by running the [regression.sh](./bin/regression.sh) bash script in the top level GSM directory, e.g.,

```bash
cd "<top level GSM directory>"
sh -e test/regression/bin/regression.sh
```

Executing this will perform all regression simulations, and in addition save results to a folder labelled [resultsCurrent]. After successful build and execution, developers will want to save these results as the base line.

> Anytime regression indicates a change in numerical precision, the regression suite should be base lined after verifying if the numerical changes are intentional or not.


### Baselining Simulation Results

Simulation results will need to be baselined anytime there is intentional regression changes. Developers should take great care to verify any regression differences obvserved.

GSM results may be baselined by executing the [updateOutputs.sh](./bin/updateOutputs.sh) script. This script will overwrite the existing result set with the latest result set.

```shell
cd "top level GSM directory>"
sh -e test/regression/bin/updateOutputs.sh
```
