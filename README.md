# hurdle_ordered
This is a function to run a hurdle ordinal model. "Don't knows" are modeled using a bernoulli distribution, while the ordinal scale is modeled as ordinal.

This is still in development, so it may or may not work well. __Use at your own risk__.

## files
There are three main files that matter:
* helper_functions_from_brms.R contains a bunch of helper functions pulled from brms. These are necessary to run the posterior predictive checks, loo, etc.  
* hurdle_ordered_function.R contains the functions needed to run the model and for post-processing.  
* test with data.R contains some simulated data and runs the model. The model works, but needs more tests

## Important!!
* The formula must be wrapped in bf().
* Priors must be manually set for `hu`. I intend to eventually set default priors, once I figure out a way to do it properly.
* `bf(y | thres(num_thres, group) ~ ...)` doesn't work. I need to figure out how to properly specify the grouping structure and update the brms model.
* Thresholds are set to be flexible. Other options don't work at the moment.
* Where possible, set the "Don't knows" to be one leve higher than the ordinal scale. This makes it easier to visualize the posterior predictive checks. If you don't care what they look like, set "Don't know" to whatever you want.
* __Use at your own risk__


## Credit where credit is due
* Check out the [Stan version](https://octavio.me/posts/2021-11-21-dealing-with-dont-knows/) and [brms version](https://discourse.mc-stan.org/t/modeling-non-response-in-ordinal-survey-data/25151/2) by Octavio Medina
* The function to fit the brms model, extract and edit the Stan code, run the cmdstanr model, and put the results back in the brmsfit object comes from @ajnafa. 

Good luck. Comments and suggestions (and PRs...) welcome.
