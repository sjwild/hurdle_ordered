


cmstanr_to_brms <- function(
  .formula, 
  .nthres = NULL, 
  .DK,
  .family, .priors, .data, .chains, .cores = parallel::detectCores(logical = F),
  .warmup, 
  #.stan_vars = NULL,
  .sampling, .seed, .base_name = "/cmdtanr_file", .base_dir = getwd(), .dest = NULL, 
  .stan_dir = str_c(.base_dir, "/stan"), .metadata = FALSE, ...
) {
  
  if(is.null(.nthres)){
    print("You must manually set the number of thresholds")
  }
  
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
  stan_vars <- stanvar(scode = stan_funs, block = "functions") #+ stanvar(scode = paste0("int DK = ", .DK), block = "tdata")
  
  
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
  
  stan_data$nthres <- .nthres
  stan_code <- stancode(brmsfit_fit)
  stan_code <- gsub("y == 99", paste0("y == ", .DK), stan_code)
  stan_code <- gsub("return 0;", paste0("return ", .DK, ";"), stan_code)
  
  
  brmsfit_fit$family$thres <- data.frame(thres = 1:.nthres,
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

# Stan code
# here comes the fun
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
      return DK;
    }
    else {
      return ordered_logistic_rng(mu, c);
    }
  }
  
"

stan_var <- stanvar(scode = stan_funs, block = "functions") #+ stanvar(scode = paste0("int DK = ", .DK), block = "tdata")

log_lik_hurdle_ordinal <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i = i)
  hu <- inv_logit(get_dpar(prep, "mu", i = i))
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  eta <- thres - mu
  y <- prep$data$Y[i]
  if (y == 1L) {
    out <- log_cdf(eta[, 1L], "logit") + 
      dbinom(0, size = 1, prob = hu, log = TRUE)
  } else if (y == nthres + 1L) {
    out <- log_ccdf(eta[, y - 1L], "logit") + 
      dbinom(0, size = 1, prob = hu, log = TRUE)
  } else if (y == nthres + 2L) {
    dbinom(1, size = 1, prob = hu, log = TRUE)
  }
  else {
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
  DK <- as.numeric(max(prep$cats))
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



out <- cmstanr_to_brms(
  .formula = bf(y ~ X1 + X2 + X3 + X4 + X5,
     hu ~ 0 + X1 + X2 + X3 + X4 + X5),
  .nthres = 2, 
  .DK = 4,
  .family = hurdle_ordinal,
  .prior = c(prior(normal(0, 2), class = b)),
  .data = df,
  .chains = 4,
  .cores = 4,
  .warmup = 1000,
  .sampling = 1000,
  .seed = 1234
)

expose_functions(out, vectorize = TRUE)


get_thres <- function(prep){
  return(prep$thres$thres)
}

first_greater <- function(A, target, i = 1) {
  ifelse(target <= A[, i] | ncol(A) == i, i, first_greater(A, target, i + 1))
}

pordinal <- function(q, eta, thres, disc = 1, family = NULL, link = "logit") {
  family <- as_one_character(family)
  link <- as_one_character(link)
  args <- nlist(x = seq_len(max(q)), eta, thres, disc, link)
  p <- do_call(paste0("d", family), args)
  .fun <- function(j) rowSums(as.matrix(p[, 1:j, drop = FALSE]))
  cblapply(q, .fun)
}

cblapply <- function(X, FUN, ...) {
  do.call(cbind, lapply(X, FUN, ...))
}

as_one_character <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.character(x)
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse_combine(s, max_char = 100L)
    stop2("Cannot coerce '", s, "' to a single character value.")
  }
  x
}

nlist <- function(...) {
  m <- match.call()
  dots <- list(...)
  no_names <- is.null(names(dots))
  has_name <- if (no_names) FALSE else nzchar(names(dots))
  if (all(has_name)) return(dots)
  nms <- as.character(m)[-1]
  if (no_names) {
    names(dots) <- nms
  } else {
    names(dots)[!has_name] <- nms[!has_name]
  }
  dots
}

subset_thres <- function(prep, i) {
  thres <- prep$thres$thres
  Jthres <- prep$thres$Jthres
  if (!is.null(Jthres)) {
    thres <- thres[, Jthres[i, 1]:Jthres[i, 2], drop = FALSE]
  }
  thres
}

dcumulative <- function(x, eta, thres, disc = 1, link = "logit") {
  eta <- disc * (thres - eta)
  if (link == "identity") {
    out <- eta
  } else {
    out <- inv_link_cumulative(eta, link = link)
  }
  out[, x, drop = FALSE]
}

inv_link_cumulative <- function(x, link) {
  x <- inv_link(x, link)
  ndim <- length(dim(x))
  dim_noncat <- dim(x)[-ndim]
  ones_arr <- array(1, dim = c(dim_noncat, 1))
  zeros_arr <- array(0, dim = c(dim_noncat, 1))
  abind::abind(x, ones_arr) - abind::abind(zeros_arr, x)
}


posterior_predict_ordinal <- function(i, prep, ...) {
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  p <- pordinal(
    seq_len(nthres + 1),
    eta = get_dpar(prep, "mu", i = i),
    disc = get_dpar(prep, "disc", i = i),
    thres = thres,
    family = prep$family$family,
    link = prep$family$link
  )
  first_greater(p, target = runif(prep$ndraws, min = 0, max = 1))
}


inv_link <- function(x, link) {
  switch(link,
         identity = x,
         log = exp(x),
         logm1 = expp1(x),
         log1p = expm1(x),
         inverse = 1 / x,
         sqrt = x^2,
         "1/mu^2" = 1 / sqrt(x),
         tan_half = 2 * atan(x),
         logit = inv_logit(x),
         probit = pnorm(x),
         cauchit = pcauchy(x),
         cloglog = inv_cloglog(x),
         probit_approx = pnorm(x),
         softplus = log1p_exp(x),
         squareplus = (x + sqrt(x^2 + 4)) / 2,
         softit = inv_softit(x),
         stop2("Link '", link, "' is not supported.")
  )
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

log_diff_exp <- function(x, y) {
  stopifnot(length(x) == length(y))
  ifelse(x > y, log(exp(x) - exp(y)), NaN)
}


log_lik_weight <- function(x, i, prep) {
  weight <- prep$data$weights[i]
  if (!is.null(weight)) {
    x <- x * weight
  }
  x
}

log_cdf <- function(x, link) {
  switch(link,
         logit = log_inv_logit(x),
         probit = pnorm(x, log.p = TRUE),
         cauchit = pcauchy(x, log.p = TRUE),
         cloglog = log1m_exp(-exp(x)),
         probit_approx = pnorm(x, log.p = TRUE),
         softit = log_inv_softit(x),
         stop2("Link '", link, "' is not supported.")
  )
}

log_ccdf <- function(x, link) {
  switch(link,
         logit = log1m_inv_logit(x),
         probit = pnorm(x, log.p = TRUE, lower.tail = FALSE),
         cauchit = pcauchy(x, log.p = TRUE, lower.tail = FALSE),
         cloglog = -exp(x),
         probit_approx = pnorm(x, log.p = TRUE, lower.tail = FALSE),
         softit = log1m_inv_softit(x),
         stop2("Link '", link, "' is not supported.")
  )
}


log_inv_logit <- function(x) {
  log(inv_logit(x))
}

log1m_inv_logit <- function(x) {
  log(1 - inv_logit(x))
}
