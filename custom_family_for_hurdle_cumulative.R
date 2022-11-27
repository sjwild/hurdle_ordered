
stan_funs <- "
  real hurdle_cumulative_lpmf(int y, real mu, real hu, vector c, int vint1) { 
  //real hurdle_cumulative_lpmf(int y, vector mu, vector hu, vector c, vint1 DK) { 

    if (y == vint1) { 
      return bernoulli_lpmf(1 | hu); 
    } else { 
      return bernoulli_lpmf(0 | hu) +  
             ordered_logistic_lpmf(y | logit(mu), c); 

    } 
  }
  
  int hurdle_cumulative_rng(real hu, real mu, vector c, int vint1) {
    real p = inv_logit(hu);
    if (bernoulli_rng(p) == 1){
      return vint1;
    }
    else {
      return ordered_logistic_rng(mu, c);
    }
  }
  
  "
  
# stan_vars for custom model
stan_vars <- stanvar(scode = stan_funs, block = "functions") 


hurdle_cumulative <- 
  # Create a custom family that is logit if y = DK, cumulative if not
  custom_family("hurdle_cumulative", 
                dpars = c("mu", "hu"),
                links = c("logit", "logit"),
                specials = "ordinal",
                type = "int",
                vars = "vint1[n]",
                threshold = "flexible")


log_lik_hurdle_cumulative <- function(i, prep) {
  mu <- logit(get_dpar(prep, "mu", i = i))
  hu <- get_dpar(prep, "hu", i = i)
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  eta <- thres - mu
  y <- prep$data$Y[i]
  DK <- max(prep$data$Y)
  if (y == 1L) {
    out <- log_cdf(eta[, 1L], prep$family$link) + 
      dbinom(0, size = 1, prob = hu, log = TRUE)
  } else if (y == nthres + 1L) {
    out <- log_ccdf(eta[, y - 1L], prep$family$link) + 
      dbinom(0, size = 1, prob = hu, log = TRUE)
  } else if (y == DK) {
    out <- dbinom(1, size = 1, prob = hu, log = TRUE)
  } else {
    out <- log_diff_exp(
      log_cdf(eta[, y], prep$family$link),
      log_cdf(eta[, y - 1L], prep$family$link) 
    ) + dbinom(0, size = 1, prob = hu, log = TRUE)
  }
  log_lik_weight(out, i = i, prep = prep)
}

posterior_predict_hurdle_cumulative <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  hu <- brms::get_dpar(prep, "hu", i = i)
  thres <- subset_thres(prep)
  nthres <- NCOL(thres)
  ndraws <- prep$ndraws
  DK <- as.numeric(max(prep$data$Y))
  p <- pordinal(
    seq_len(nthres + 1),
    eta = logit(mu),
    disc = 1,
    thres = thres,
    family = "cumulative",
    link = prep$family$link
  )
  draws <- first_greater(p, target = runif(prep$ndraws, min = 0, max = 1))
  theta <- runif(ndraws, 0, 1)
  draws[hu > theta] <- DK
  return(draws)
}


posterior_epred_hurdle_cumulative <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  hu <- brms::get_dpar(prep, "hu")
  return(mu * hu)
}


