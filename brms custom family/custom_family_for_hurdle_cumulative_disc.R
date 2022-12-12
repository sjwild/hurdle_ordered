
hurdle_cumulative <- 
  # Create a custom family that is logit if y = DK, cumulative if not
  custom_family("hurdle_cumulative", 
                dpars = c("mu", "hu", "disc"),
                links = c("logit", "logit", "log"),
                specials = "ordinal",
                type = "int",
                vars = "vint1[n]",
                threshold = "flexible")


stan_funs <- "
  real hurdle_cumulative_lpmf(int y, real mu, real hu, real disc, vector thres, int vint1) { 
  
  int nthres = num_elements(c);
  real mu_logit = logit(mu);
  
    if (y == vint1) { 
      return bernoulli_lpmf(1 | hu); 
    } else { 
      if (y == 1) {
        return(log_inv_logit(disc * (thres[1] - mu_logit)) +
        bernoulli_lpmf(0 | hu);
      } else if (y == nthres + 1) {
        return log1m_inv_logit(disc * (thres[nthres] - mu_logit)) +
        bernoulli_lpmf(0 | hu);
      } else {
         return log_diff_exp(
           log_inv_logit(disc * (thres[y] - mu_logit)), 
           log_inv_logit(disc * (thres[y - 1] - mu_logit))
       ) + bernoulli_lpmf(0 | hu);
      }
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

log_lik_hurdle_cumulative <- function(i, prep) {
  mu <- logit(get_dpar(prep, "mu", i = i))
  hu <- get_dpar(prep, "hu", i = i)
  disc <- get_dpar(prep, "disc", i = i)
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  eta <- disc * (thres - mu)
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
  disc <- brms::get_dpar(prep, "disc", i = i)
  thres <- subset_thres(prep)
  nthres <- NCOL(thres)
  ndraws <- prep$ndraws
  DK <- prep$data$vint1
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
  draws[theta <- hu] <- DK[1]
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



#### function to run brms model ####
cmstanr_to_brms <- function(
  .formula,  # must be in the form of bf(y ~ ., hu ~ .)
  .outcome = NULL, # put the outcome here. It is necessary to get the number of thresholds
  .family = hurdle_cumulative,
  .priors, # if using default priors, brms will crash. You must adjust the prior for hu, as brms use beta as the default
  .data, # dataframe here
  .chains = 4, # 4 chains by default
  .cores = parallel::detectCores(logical = F),
  .warmup = 1000, # 1000 warmup iterations
  .sampling = 1000, # 1000 sampling iterations
  .seed, 
  .base_name = "/cmdtanr_file", 
  .base_dir = here(), 
  .dest = NULL, 
  .stan_dir = str_c(.base_dir, "/stan"), 
  .metadata = FALSE, 
  ...
) {
  
  if(is.null(.outcome)){
    stop("You must input the outcome. Do not forget to set DK to be one higher than the scale")
  }
  
  nthres <- length(unique(.outcome)) - 2
  DK <- max(.outcome)
  
  # stan funs for custom family
  stan_funs <- "
  real hurdle_cumulative_lpmf(int y, real mu, real hu, real disc, vector thres, int vint1) { 
  
  int nthres = num_elements(thres);
  real mu_logit = logit(mu);
  
    if (y == vint1) { 
      return bernoulli_lpmf(1 | hu); 
    } else { 
      if (y == 1) {
        return log_inv_logit(disc * (thres[1] - mu_logit)) +
        bernoulli_lpmf(0 | hu);
      } else if (y == nthres + 1) {
        return log1m_inv_logit(disc * (thres[nthres] - mu_logit)) +
        bernoulli_lpmf(0 | hu);
      } else {
         return log_diff_exp(
           log_inv_logit(disc * (thres[y] - mu_logit)), 
           log_inv_logit(disc * (thres[y - 1] - mu_logit))
       ) + bernoulli_lpmf(0 | hu);
      }
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
  
  
  # Generate a temporary brmsfit object
  brmsfit_fit <- brm(
    formula = .formula, # Model formula for the brmsfit object
    family = .family, # Model likelihood
    prior = .priors, # Priors for the brms model
    data = .data, # Data for the model parameters
    chains = 0, # Number of MCMC Chains
    cores = .cores, # Number of Physical Cores
    seed = .seed, # Random Number Seed for Reproducibility
    backend = "cmdstanr", # Requires cmdstanr as a backend
    empty = TRUE, # Do not fit the model object
    save_pars = save_pars(all = TRUE), # Save all model parameters
    stanvars = stan_vars
  )
  
  # Generate the stan data for the model
  stan_data <- make_standata(
    formula = .formula, # Model formula from the brmsfit object
    family = .family, # Model likelihood
    prior = .priors, # Priors from the brms model
    data = .data, # Data for the model parameters
    stanvars = stan_vars
  )
  
  stan_data$nthres <- nthres
  stan_code <- stancode(brmsfit_fit)
  
  
  brmsfit_fit$family$thres <- data.frame(thres = 1:nthres,
                                         #group = as.character(c(1:(.nthres))),
                                         group = ""
  )
  
  
  # Generate the cmdstan model
  cmdstanr_mod <- cmdstan_model(
    stan_file = write_stan_file(stan_code),
    dir = str_c(.base_dir, .dest),
    force_recompile = TRUE
  )
  
  # Fit the cmdstan model
  cmdstanr_mod_fit <- cmdstanr_mod$sample(
    data = stan_data, # A list object to pass the data to Stan
    seed = .seed, # Random Number Seed for Reproducibility
    output_dir = str_c(.base_dir, .dest), # Output directory to write the csv files
    output_basename = .base_name, # Prefix for the csv files
    parallel_chains = .chains, # Number of chains to run in parallel
    chains = .chains, # Number of chains to run
    iter_warmup = .warmup, # Warmup Iterations
    iter_sampling = .sampling, # Sampling Iterations
    ... # Additional arguments for the sampler
  )
  
  # Save the cmdstanr environment
  cmdstanr_mod_fit$save_object(file = str_c(.base_dir, .dest, .base_name, ".rds"))
  
  # Convert the environment to a stanfit object
  stanfit_mod <- rstan::read_stan_csv(cmdstanr_mod_fit$output_files())
  
  # Add the stanfit object to the brmsfit object
  brmsfit_fit$fit <- stanfit_mod
  
  # Rename the parameters
  brmsfit_fit <- rename_pars(brmsfit_fit)
  
  # add in the stan code
  brmsfit_fit$model <- stan_code
  
  # Add any Profiles
  #brmsfit_fit$profiles <- bind_rows(cmdstanr_mod_fit$profiles(), .id = "chain")
  
  # Construct a field with relevant metadata if .metadata is TRUE
  if (.metadata == TRUE){
    brmsfit_fit$meta_data <- map_dfr(
      .x = brmsfit_fit$fit@stan_args,
      .f = ~ tibble(
        warmup = str_remove_all(.x$time_info[1], "[^[:digit:]|\\.]"),
        sampling = str_remove_all(.x$time_info[2], "[^[:digit:]|\\.]"),
        total = str_remove_all(.x$time_info[3], "[^[:digit:]|\\.]"),
        misc = str_remove_all(.x$time_info[4], "[^[:digit:]|\\.]"),
        metadata = c(
          str_c("stanc_version:", .x$stanc_version[1], sep = " "),
          str_c("opencl_device_name:", .x$opencl_device_name[1], sep = " "),
          str_c("opencl_platform_name", .x$opencl_platform_name[1], sep = " "),
          str_c("date", .x$start_datetime[1], sep = " ")
        )
      ),
      .id = "chain"
    )
  }
  
  # Write the brmsfit object to an rds file
  write_rds(brmsfit_fit, str_c(.base_dir, .dest, .base_name, "_brms.rds"))
  
  # Return the brms object
  return(brmsfit_fit)
}



