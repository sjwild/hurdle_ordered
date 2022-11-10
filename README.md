# hurdle_ordered
This is a function to run a hurdle ordinal model. "Don't knows" are modeled using a bernoulli distribution, while the ordinal scale is modeled as ordinal.

This is still in development, so it may or may not work well. Use at your own risk.

## files
There are three files:
* helper_functions_from_brms.R contains a bunch of helper functions pulled from brms. These are necessary to run the posterior predictive checks, loo, etc.  
* hurdle_ordered_function.R contains the functions needed to run the model and for post-processing.  
* test with data.R contains some simulated data and runs the model. The model works, but needs more tests

## Important!!
* The formula must be wrapped in bf().
* Priors must be manually set. I intend to eventually set default priors

Good luck.
