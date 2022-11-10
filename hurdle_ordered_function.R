


cmstanr_to_brms <- function(
  .formula, 
  .outcome = NULL, 
  #.DK,
  .family, .priors, .data, .chains, .cores = parallel::detectCores(logical = F),
  .warmup, 
  #.stan_vars = NULL,
  .sampling, .seed, .base_name = "/cmdtanr_file", .base_dir = here(), .dest = NULL, 
  .stan_dir = str_c(.base_dir, "/stan"), .metadata = FALSE, ...
) {
  
  if(is.null(.outcome)){
    stop("You must input the outcome. Do not forget to set DK to be one higher than the scale")
  }
  
  nthres <- length(unique(.outcome)) - 2
  DK <- max(.outcome)
  
  # stan funs for custom family
  stan_funs <- "
  real hurdle_ordinal_lpmf(int y, real mu, real hu, vector c) { 
  //real hurdle_ordinal_lpmf(int y, vector mu, vector hu, vector c, vint1 DK) { 

    if (y == 99) { 
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
  stan_code <- gsub("y == 99", paste0("y == ", DK), stan_code)
  stan_code <- gsub("return 0;", paste0("return ", DK, ";"), stan_code)
  
  
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


hurdle_ordinal <- 
  # Create a custom family that is logit if y = 0, normal/gaussian if not
  custom_family("hurdle_ordinal", 
                dpars = c("mu", "hu"),
                links = c("identity", "identity"),
                specials = "ordinal",
                type = "int",
                threshold = "flexible")

log_lik_hurdle_ordinal <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i = i)
  hu <- inv_logit(get_dpar(prep, "hu", i = i))
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  eta <- thres - mu
  y <- prep$data$Y[i]
  DK <- max(prep$data$Y)
  if (y == 1L) {
    out <- log_cdf(eta[, 1L], "logit") + 
      dbinom(0, size = 1, prob = hu, log = TRUE)
  } else if (y == nthres + 1L) {
    out <- log_ccdf(eta[, y - 1L], "logit") + 
      dbinom(0, size = 1, prob = hu, log = TRUE)
  } else if (y == DK) {
    out <- dbinom(1, size = 1, prob = hu, log = TRUE)
  } else {
    out <- log_diff_exp(
      log_cdf(eta[, y], "logit"),
      log_cdf(eta[, y - 1L], "logit") 
    ) + dbinom(0, size = 1, prob = hu, log = TRUE)
  }
  log_lik_weight(out, i = i, prep = prep)
}

posterior_predict_hurdle_ordinal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  hu <- brms::get_dpar(prep, "hu", i = i)
  thres <- subset_thres(prep)
  nthres <- NCOL(thres)
  ndraws <- prep$ndraws
  DK <- as.numeric(max(prep$data$Y))
  p <- pordinal(
    seq_len(nthres + 1),
    eta = get_dpar(prep, "mu", i = i),
    disc = 1,
    thres = thres,
    family = "cumulative",
    link = "logit"
  )
  draws <- first_greater(p, target = runif(prep$ndraws, min = 0, max = 1))
  theta <- runif(ndraws, 0, 1)
  draws[inv_logit(hu) > theta] <- DK
  return(draws)
}


posterior_epred_hurdle_ordinal <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  hu <- brms::get_dpar(prep, "hu")
  return(inv_logit(mu) * inv_logit(hu))
}


