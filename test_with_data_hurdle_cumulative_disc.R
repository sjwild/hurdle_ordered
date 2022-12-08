library(brms)
library(cmdstanr)
library(posterior)
library(magrittr)
library(tidyverse)
library(here)

source("helper_functions_from_brms.R")
source("custom_family_for_hurdle_cumulative_disc.R")


#### Simulate some data for testing ####
N <- 10000
P <- 5
cutpoints <- c(1, -1)
y <- rep(0, N)


set.seed(143140914)
X <- matrix(runif(N * P), N, P)
betas <- c(-.15, -.5, 3, .5, -1)
gammas <- c(3, -2, -2, -1, -.5)
p2 <- plogis(cutpoints[1] + X %*% betas)
p3 <- plogis(cutpoints[2] + X %*% betas)
p1 <- 1 - p2
p2 <- p2 - p3
for(i in 1:N){
  y[i] <- sample(1:3, 
                 size = 1,
                 prob = c(p1[i], p2[i], p3[i]))
}


# Simulate non-response bias
p <- plogis(X %*% gammas)
dk <- rbinom(N, 1, p)
y[dk == 1] = 99
summary(as.factor(y))

df <- data.frame(y = y, X)



# test function
out <- cmstanr_to_brms(
  .formula = bf(y | vint(99) ~ X1 + X2 + X3 + X4 + X5,
                hu ~ 0 + X1 + X2 + X3 + X4 + X5),
  .outcome = df$y,
  .family = hurdle_cumulative,
  .prior = c(prior(normal(0, 2), class = b)),
  .data = df,
  .chains = 4,
  .cores = 4,
  .warmup = 1000,
  .sampling = 1000,
  .seed = 1234
)

expose_functions(out, vectorize = TRUE)


# model checks
summary(out)
pp_check(out, type = "bars") + theme_minimal()
loo(out)



merkel <- read.csv("https://raw.githubusercontent.com/octmedina/zi-ordinal/main/merkel_data.csv")
merkel$confid_merkel[merkel$confid_merkel == 0] <- 97

out_merkel_disc <- cmstanr_to_brms(
  .formula = bf(confid_merkel | vint(97) ~ edu + race + income + party,
                #disc ~ 1 + edu + race + income + party,
                hu ~ 1 + edu + race + income + party
  ),
  .outcome = merkel$confid_merkel,
  .family = hurdle_cumulative2,
  .prior = c(prior(normal(0, 2), class = b),
             prior(normal(0, 2), class = Intercept, dpar = hu)),
  .data = merkel,
  .chains = 4,
  .cores = 4,
  .warmup = 1000,
  .sampling = 2000,
  #adapt_delta = 0.99, # if testing disc in formula
  .seed = 1234
)

summary(out_merkel_disc)
pp_check(out_merkel_disc, type = "bars") + theme_minimal()
loo(out_merkel_disc)