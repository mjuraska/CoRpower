# CoRpower 1.0.5

* Removed the dependency on the package osDesign, which is scheduled for archival on CRAN. The function `tps()` from osDesign (version 1.8) is now included as an internal function with unchanged numerical code.
* The license changed to GPL (>= 3) to match that of the included osDesign code.
* `computePower()` with `tpsMethod = "WL"` no longer fails: the WL fit provides only the empirical variance estimate, which is now used for the Wald test.
* Fixed handling of a dichotomous biomarker that caused crashing.
* Fixed invalid allocation of simulated cases among latent subgroups in the trichotomous biomarker simulation.
