# Supplementary
This code reproduces all results in the corncob manuscript.

Simulations were conducted using the (highly recommended!) [simulator R package](http://faculty.bscb.cornell.edu/~bien/simulator.html) (Bien, 2016).

Parts of the parallelization are specific to the computing systems used.
To run the code on any computer:
 * For Figures 2 and 3 - remove the parallel component of calls to `run_method()` in main.R.
 * For Figures 4 and 5 - remove any calls dependent on the parallel package, and switch `parLapply` to `lapply`.
