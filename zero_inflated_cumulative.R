
stan_funs <- "
  real zi_cumulative_lpmf(int y, real mu, real zi, vector c, int vint1) { 

    if (y == vint1) { 
      return log_sum_exp(bernoulli_lpmf(1 | zi),  
                         bernoulli_lpmf(0 | zi) + 
                         ordered_logistic_lpmf(y | logit(mu), c)); 
    } else { 
      return bernoulli_lpmf(0 | zi) +  
             ordered_logistic_lpmf(y | logit(mu), c); 
    } 
  }
  
  "

zi_cumulative <- 
  # Create a custom family that is logit if y = DK, cumulative if not
  custom_family("zi_cumulative", 
                dpars = c("mu", "zi"),
                links = c("logit", "logit"),
                specials = "ordinal",
                type = "int",
                vars = "vint1[n]",
                threshold = "flexible")


log_lik_zi_cumulative <- function(i, prep) {
  mu <- logit(get_dpar(prep, "mu", i = i))
  theta <- get_dpar(prep, "zi", i = i)
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  eta <- thres - mu
  y <- prep$data$Y[i]
  zi <- prep$data$vint1[i]
  if (y == 1L & y == zi) {
    out <- log_cdf(eta[, 1L], prep$family$link) + 
      dbinom(1, size = 1, prob = theta, log = TRUE)
  } else if (y == 1L) {
    out <- log_cdf(eta[, 1L], prep$family$link) + 
      dbinom(0, size = 1, prob = theta, log = TRUE)
  } else if (y == nthres + 1L & y == zi) {
    out <- log_ccdf(eta[, y - 1L], prep$family$link) + 
      dbinom(1, size = 1, prob = theta, log = TRUE)
  } else if (y == nthres + 1L) {
    out <- log_ccdf(eta[, y - 1L], prep$family$link) + 
      dbinom(0, size = 1, prob = theta, log = TRUE)
  } else if (y == zi) {
    out <- log_diff_exp(
      log_cdf(eta[, y], prep$family$link),
      log_cdf(eta[, y - 1L], prep$family$link) 
    ) + dbinom(1, size = 1, prob = theta, log = TRUE)
  } else {
    out <- log_diff_exp(
      log_cdf(eta[, y], prep$family$link),
      log_cdf(eta[, y - 1L], prep$family$link) 
    ) + dbinom(0, size = 1, prob = theta, log = TRUE)
  }
  log_lik_weight(out, i = i, prep = prep)
}




# Gotta work on this
posterior_predict_zi_cumulative <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  theta <- brms::get_dpar(prep, "zi", i = i)
  thres <- subset_thres(prep)
  nthres <- NCOL(thres)
  ndraws <- prep$ndraws
  ZI <- prep$data$vint1
  p <- pordinal(
    seq_len(nthres + 1),
    eta = logit(mu),
    disc = 1,
    thres = thres,
    family = "cumulative",
    link = "logit"
  )
  draws <- first_greater(p, target = runif(prep$ndraws, min = 0, max = 1))
  theta_draws <- runif(ndraws, 0, 1)
  draws[theta_draws < theta] <- ZI[1]
  return(draws)
}

#and gotta work on this
posterior_epred_zi_cumulative <- function(prep) {
  #dens <- get(paste0("d", "cumulative"), mode = "function")
  # the linear scale has one column less than the response scale
  adjust <- ifelse(prep$family$link == "identity", 0, 1)
  #ncat_max <- max(prep$data$nthres) + adjust
  ncat_max <- NCOL(prep$thres$thres) + adjust
  #nact_min <- min(prep$data$nthres) + adjust
  #init_mat <- matrix(ifelse(prep$family$link == "identity", NA, 0),
  #                   nrow = prep$ndraws,
  #                   ncol = ncat_max - nact_min)
  args <- list(link = prep$family$link)
  out <- vector("list", prep$nobs)
  ZI <- prep$data$vint1[1]
  
  #for (i in seq_along(out)) {
  for (i in seq_along(out)) {
    args_i <- args
    args_i$eta <- logit(slice_col(get_dpar(prep, "mu", i)))
    #args_i$disc <- slice_col(get_dpar(prep, "disc", i))
    args_i$disc <- 1
    args_i$thres <- subset_thres(prep, i)
    ncat_i <- NCOL(args_i$thres) + adjust
    args_i$x <- seq_len(ncat_i)
    out[[i]] <- do_call(dcumulative, args_i)
    zi <- get_dpar(prep, "zi", i)
    out[[i]] <- out[[i]] * (1 - zi)
    out[[i]][,ZI] <- out[[i]][,ZI] + zi
  
    # this section not needed right now because we cannot use thresh
    # But can leave as is because 
    #if (ncat_i < ncat_max) {
    #  sel <- seq_len(ncat_max - ncat_i)
    #  out[[i]] <- cbind(out[[i]], init_mat[, sel])
    #}
  }
  
  out <- abind(out, along = 3)
  out <- aperm(out, perm = c(1, 3, 2))
  dimnames(out)[[3]] <- seq_len(ncat_max)
  return(out)
  
}



# test with data
merkel <- read.csv("https://raw.githubusercontent.com/octmedina/zi-ordinal/main/merkel_data.csv")
merkel$confid_merkel[merkel$confid_merkel == 0] <- 3

out_merkel_zi <- brm(bf(confid_merkel | vint(3) ~ edu + race + income + party,
                #disc ~ 1 + edu,
                zi ~ 1 + edu + race + income + party),
                data = merkel,
                family = zi_cumulative,
                stan_funs = stan_funs,
                cores = 4,
                chains = 4,
                iter = 4000,
                warmup = 2000,
                backend = "cmdstanr")
                
loo(out_merkel_zi)
pp_check(out_merkel_zi, type = "bars", ndraws = 50)        
summary(out_merkel_zi)





#### Simulate some data for testing ####
N <- 10000
N_zi <- 1000
P <- 5
cutpoints <- c(1, -1)
y_zi <- rep(0, N)


set.seed(143140914)
X <- matrix(runif(N * P), N, P)
betas <- c(-.15, -.5, 3, .5, -1)
gammas <- c(3, -2, -2, -1, -.5)
p2 <- plogis(cutpoints[1] + X %*% betas)
p3 <- plogis(cutpoints[2] + X %*% betas)
p1 <- 1 - p2
p2 <- p2 - p3
for(i in 1:N){
  y_zi[i] <- sample(1:3, 
                 size = 1,
                 prob = c(p1[i], p2[i], p3[i]))
}


# Simulate non-response bias
p <- plogis(X %*% gammas)
dk <- rbinom(N, 1, p)
y_zi[dk == 1] = 2

df_zi <- data.frame(y = y_zi, X)


out_sim_zi <- brm(bf(y | vint(2) ~ X1 + X2 + X3 + X4 + X5,
                  zi ~ 0 + X1 + X2 + X3 + X4 + X5),
                  prior = c(prior(normal(0, 2), class = b)),
                  data = df_zi,
                  family = zi_cumulative,
                  stan_funs = stan_funs,
                  cores = 4,
                  chains = 4,
                  iter = 2000,
                  warmup = 1000,
                  backend = "cmdstanr")
summary(out_sim_zi)                  
pp_check(out_sim_zi, type = "bars")

