
stan_funs <- "
  real hurdle_cumulative_lpmf(int y, real mu, real hu, vector c, int vint1) { 

    if (y == vint1) { 
      return bernoulli_lpmf(1 | hu); 
    } else { 
      return bernoulli_lpmf(0 | hu) +  
             ordered_logistic_lpmf(y | logit(mu), c); 
    } 
  }
  
  int hurdle_cumulative_rng(real mu, real hu, vector c, int vint1) {
    if (bernoulli_rng(hu) == 1){
      return vint1;
    }
    else {
      return ordered_logistic_rng(logit(mu), c);
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
  DK <- prep$data$vint1[i]
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
  #disc <- brms::get_dpar(prep, "disc", i = i)
  disc = 1
  thres <- subset_thres(prep)
  nthres <- NCOL(thres)
  ndraws <- prep$ndraws
  DK <- as.numeric(max(prep$data$Y))
  p <- pordinal(
    seq_len(nthres + 1),
    eta = logit(mu),
    disc = disc,
    thres = thres,
    family = "cumulative",
    link = prep$family$link
  )
  draws <- first_greater(p, target = runif(prep$ndraws, min = 0, max = 1))
  theta <- runif(ndraws, 0, 1)
  draws[hu > theta] <- DK[draws[hu > theta]]
  return(draws)
}



posterior_epred_hurdle_cumulative <- function(prep) {
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
  
  #for (i in seq_along(out)) {
  for (i in seq_along(out)) {
    args_i <- args
    args_i$eta <- logit(slice_col(get_dpar(prep, "mu", i)))
    args_i$disc <- slice_col(get_dpar(prep, "disc", i))
    args_i$thres <- subset_thres(prep, i)
    ncat_i <- NCOL(args_i$thres) + adjust
    args_i$x <- seq_len(ncat_i)
    out[[i]] <- do_call(dcumulative, args_i)
    
    # this section not needed right now because we cannot use thresh
    # But can leave as is because 
    #if (ncat_i < ncat_max) {
    #  sel <- seq_len(ncat_max - ncat_i)
    #  out[[i]] <- cbind(out[[i]], init_mat[, sel])
    #}
    
    hu <- get_dpar(prep, "hu", i)
    out[[i]] <- cbind(out[[i]] * (1 - hu), hu)
    
  }
  
  out <- abind(out, along = 3)
  out <- aperm(out, perm = c(1, 3, 2))
  DK <- prep$data$vint1
  dimnames(out)[[3]] <- c(seq_len(ncat_max), paste0(DK[1]))
  return(out)
  
}





