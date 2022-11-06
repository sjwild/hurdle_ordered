library(brms)
library(cmdstanr)
library(posterior)
library(magrittr)
library(tidyverse)

# simulate simple data
N <- 1000
prob_dk <- .1
probs = c(0.2, 0.1, 0.35, 0.2, 0.15)

set.seed(984301483)
dk <- rbinom(N, 1, prob_dk)
y <- sample(1:5, 
            N, 
            replace = TRUE, 
            prob = probs) * (1 - dk)
y[y == 0] = 99 # DK usually coded as 99



# run model
hu_ordered <- cmdstan_model("hurdle_ordered.stan")

stan_dat <- list(
  N_cutpoints = 5,
  N_obs = N,
  P = 1,
  y = y
  
)

fit <- hu_ordered$sample(data = stan_dat,
                         seed = 134123,
                         chains = 4,
                         parallel_chains = 4,
                         refresh = 100)

fit$summary()



# more complex version
N <- 10000
P <- 5
cutpoints <- c(1, -1)
y_star <- y <- rep(0, N)


set.seed(143140914)
X <- matrix(runif(N * P), N, P)
betas <- c(-.15, -.5, 3, .5, -1)
gammas <- c(3, -2, -2, -1, -.5)
p2 <- plogis(cutpoints[1] + X %*% betas)
p3 <- plogis(cutpoints[2] + X %*% betas)
p1 <- 1 - p2
p2 <- p2 - p3
for(i in 1:N){
  y_star[i] <- y[i] <- sample(1:3, 
                              size = 1,
                              prob = c(p1[i], p2[i], p3[i]))
}


# Simulate non-response bias
p <- plogis(X %*% gammas)
dk <- rbinom(N, 1, p)
y[dk == 1] = 4
summary(as.factor(y))

# new stan
hu_ordered_2 <- cmdstan_model("hurdle_ordered_2.stan")

stan_dat_2 <- list(
  N_levels = 3,
  N_obs = N,
  P = 5,
  y = y,
  X = X,
  DK = 4
  
)

fit_2 <- hu_ordered_2$sample(data = stan_dat_2,
                         seed = 1323,
                         chains = 4,
                         parallel_chains = 4,
                         refresh = 100)

fit_2$summary(c("beta", "cutpoints"))


# Put in dataframe for brms
df <- data.frame(y = y, X)

mod_drop_dk <- brm(y ~ 1 + X1 + X2 + X3 + X4 + X5,
                   family = cumulative,
                   data = df[df$y != 4,],
                   prior = c(prior(normal(0, 1.5), class = b)),
                   cores = 4,
                   chains = 4,
                   backend = "cmdstanr")


fit_2$summary(c("beta", "cutpoints"))
summary(mod_drop_dk)














#### brms model ####
# custom family
hurdle_ordinal <- 
  # Create a custom family that is logit if y = 0, normal/gaussian if not
  custom_family("hurdle_ordinal", 
                dpars = c("mu", "hu"),
                links = c("identity", "identity"),
                specials = "ordinal",
                type = "int",
                threshold = "flexible")

# Stan code
# here comes the fun
stan_funs <- "
  real hurdle_ordinal_lpmf(int y, real mu, real hu, vector c) { 
  //real hurdle_ordinal_lpmf(int y, vector mu, vector hu, vector c, vint1 DK) { 

    if (y == 4) { 
      return bernoulli_logit_lpmf(1 | hu); 
    } else { 
      return bernoulli_logit_lpmf(0 | hu) +  
             ordered_logistic_lpmf(y | mu, c); 
    } 
  }
  
  int hurdle_ordinal_rng(real hu, real mu, vector c) {
    real p = inv_logit(hu);
    if (bernoulli_rng(p) == 1){
      return 0;
    }
    else {
      return ordered_logistic_rng(mu, c);
    }
  }
  
"

stan_var <- stanvar(scode = stan_funs, block = "functions")# + stanvar(scode = "int DK;", block = "parameters")


fit_3 <- brm(bf(y ~ X1 + X2 + X3 + X4 + X5,
                hu ~ 0 + X1 + X2 + X3 + X4 + X5),
             family = hurdle_ordinal,
             stanvars = stan_var,
             data = df,
             prior = c(prior(normal(0, 2), class = b)),
             chains = 0,
             cores = 4)

sc <- stancode(fit_3)

sc <- gsub("vector\\[nthres\\] Intercept", "vector\\[nthres-1\\] Intercept", sc)

fit_3$model <- sc #stan_model(model_code = sc)
fit_3$fit@par_dims$Intercept <- 2
fit_3$family$thres <- 1:2
fit_3$family$cats <- as.character(c(1:3))


fit_3 <- update(fit_3,
                cores = 4,
                chains = 4)

summary(fit_3)



brm_hurdle_ordinal <- function(formula,
                               nthres,
                               dk_value = 99,
                               data) {
  
  
  
  
  
  
  
  
}



# this feels sort of hacky, but to define the ordered c_int variable i had to
# add this to the parameters block, as well as the number of categories (n_thresh)
#ordered_var <- stanvar(scode = "ordered[n_thresh] c_int;", block = "parameters")
#ncat_var <- stanvar(x = 3, name = "n_thresh", scode = "int n_thresh;")
#stanvars <- ordered_var + ncat_var + stan_funs


# 
stan_data <- "nthres = 3
   "


# Prepare Stan code for use in brm()
stanvars <- stanvar(scode = stan_funs, block = "functions")

make_stancode(y ~ 1, family = cumulative(), data = data.frame(y = y[y != 99]))




# empty brms









# Test with Octavio data
data <- read.csv("https://raw.githubusercontent.com/octmedina/zi-ordinal/main/merkel_data.csv") 
data$confid_merkel[data$confid_merkel == 0] = 99
stan_dat_3 <- list(
  N_cutpoints = 4,
  N_obs = nrow(data),
  P = 4,
  y = data$confid_merkel,
  X = as.matrix(data[,2:5])
  
)

fit_3 <- hu_ordered_2$sample(data = stan_dat_3,
                             seed = 1323,
                             chains = 4,
                             parallel_chains = 4,
                             refresh = 100)
