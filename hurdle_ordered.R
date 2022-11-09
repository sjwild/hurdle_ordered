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
hurdle_ordinal2 <- 
  # Create a custom family that is logit if y = 0, normal/gaussian if not
  custom_family("hurdle_ordinal2", 
                dpars = c("mu", "hu", "c"),
                links = c("identity", "identity", "identity"),
                specials = "ordinal",
                type = "int",
                threshold = "flexible")

# Stan code
# here comes the fun
stan_funs <- "
  real hurdle_ordinal2_lpmf(int y, real mu, real hu, vector c, vector Intercept) { 
  //real hurdle_ordinal2_lpmf(int y, vector mu, vector hu, vector disc) { 

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

stan_tp <- "
    vector[nthres] c = Intercept;"

stan_var <- stanvar(scode = stan_funs, block = "functions") + stanvar(scode = stan_tp, block = "tparameters")

make_stancode(bf(y ~ X1 + X2 + X3 + X4 + X5,
                 hu ~ 0 + X1 + X2 + X3 + X4 + X5),
              family = hurdle_ordinal2,
              stanvars = stan_var,
              data = df,
              prior = c(prior(normal(0, 2), class = b)))








fit_3 <- brm(bf(y ~ X1 + X2 + X3 + X4 + X5,
                hu ~ 0 + X1 + X2 + X3 + X4 + X5),
             family = hurdle_ordinal2,
             stanvars = stan_var,
             data = df,
             prior = c(prior(normal(0, 2), class = b)),
             seed = 1244,
             chains = 0,
             cores = 4,
             backend = "cmdstanr")

sc <- stancode(fit_3)

sc <- gsub("vector\\[nthres\\]", "vector\\[nthres-1\\]", sc)

cmd_mod <- cmdstan_model(write_stan_file(sc))
md <- read_cmdstan_csv(cmd_mod$output_files())


cmdstan_fit <- cmd_mod$sample(data = standata(fit_3),
                              parallel_chains = 4,
                              refresh = 100,
                              seed = 1244)



metadata <- read_cmdstan_csv(cmdstan_fit$output_files(),
                             variables = "", 
                             sampler_diagnostics = "")
variables <- repair_variable_names(metadata$metadata$variables)
variables <- unique(sub("\\[.+", "", variables))
#variables <- setdiff(variables, exclude)
out <- read_csv_as_stanfit(cmdstan_fit$output_files(), variables = variables)
out <- repair_stanfit(out)
# allow updating the model without recompilation
attributes(out)$CmdStanModel <- model
attributes(out)$metadata <- metadata
if (empty_model) {
  # allow correct updating of an 'empty' model
  out@sim <- list()
}
out
}




fit_3$model <- sc 




fit_3$fit <- brms:::.fit_model_cmdstanr(model_code = sc, 
                                        sdata = standata(fit_3),
                                        algorithm = fit_3$algorithm,
                                        iter = 2000,
                                        warmup = 1000,
                                        thin = ,
                                        chains = ,
                                        cores = ,
                                        threads = fit_3$threads$threads,
                                        opencl = fit_3$opencl$ids,
                                        init = fit_3$fit@inits,
                                        future = FALSE,
                                        control = 
                                        seed = 1244)


metadata <- cmdstanr::read_cmdstan_csv(
  out$output_files(), variables = "", sampler_diagnostics = ""
)
# ensure that only relevant variables are read from CSV
variables <- repair_variable_names(metadata$metadata$variables)
variables <- unique(sub("\\[.+", "", variables))
variables <- setdiff(variables, exclude)
# transform into stanfit object for consistent output structure
out <- read_csv_as_stanfit(out$output_files(), variables = variables)
out <- repair_stanfit(out)
# allow updating the model without recompilation
attributes(out)$CmdStanModel <- model
attributes(out)$metadata <- metadata
if (empty_model) {
  # allow correct updating of an 'empty' model
  out@sim <- list()
}
out
function(model, sdata, algorithm, iter, warmup, thin,
         chains, cores, threads, opencl, init, exclude,
         seed, control, silent, future, ...)
  

fit_3$fit@par_dims$Intercept <- 2
fit_3$family$thres$thres <- 1:2
fit_3$family$cats <- as.character(c(1:3))


fit_3 <- update(fit_3,
                recompile = FALSE)

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

# functions
repair_variable_names <- function(x) {
  x <- sub("\\.", "[", x)
  x <- gsub("\\.", ",", x)
  x[grep("\\[", x)] <- paste0(x[grep("\\[", x)], "]")
  x
}

repair_stanfit <- function(x) {
  stopifnot(is.stanfit(x))
  if (!length(x@sim$fnames_oi)) {
    # nothing to rename
    return(x)
  }
  # the posterior package cannot deal with non-unique parameter names
  # this case happens rarely but might happen when sample_prior = "yes"
  x@sim$fnames_oi <- make.unique(as.character(x@sim$fnames_oi), "__")
  for (i in seq_along(x@sim$samples)) {
    # stanfit may have renamed dimension suffixes (#1218)
    if (length(x@sim$samples[[i]]) == length(x@sim$fnames_oi)) {
      names(x@sim$samples[[i]]) <- x@sim$fnames_oi
    }
  }
  x
}


read_csv_as_stanfit <- function(files, variables = NULL,
                                sampler_diagnostics = NULL) {
  
  csfit <- cmdstanr::read_cmdstan_csv(
    files = files, variables = variables,
    sampler_diagnostics = sampler_diagnostics,
    format = NULL
  )
  
  # @model_name
  model_name = gsub(".csv", "", basename(files[[1]]))
  
  # @model_pars
  svars <- csfit$metadata$stan_variables
  if (!is.null(variables)) {
    variables_main <- unique(gsub("\\[.*\\]", "", variables))
    svars <- intersect(variables_main, svars)
  }
  if ("lp__" %in% svars) {
    svars <- c(setdiff(svars, "lp__"), "lp__")
  }
  pars_oi <- svars
  par_names <- csfit$metadata$model_params
  
  # @par_dims
  par_dims <- vector("list", length(svars))
  
  names(par_dims) <- svars
  par_dims <- lapply(par_dims, function(x) x <- integer(0))
  
  pdims_num <- ulapply(
    svars, function(x) sum(grepl(paste0("^", x, "\\[.*\\]$"), par_names))
  )
  par_dims[pdims_num != 0] <-
    csfit$metadata$stan_variable_sizes[svars][pdims_num != 0]
  
  # @mode
  mode <- 0L
  
  # @sim
  rstan_diagn_order <- c("accept_stat__", "treedepth__", "stepsize__",
                         "divergent__", "n_leapfrog__", "energy__")
  
  if (!is.null(sampler_diagnostics)) {
    rstan_diagn_order <- rstan_diagn_order[rstan_diagn_order %in% sampler_diagnostics]
  }
  
  res_vars <- c(".chain", ".iteration", ".draw")
  if ("post_warmup_draws" %in% names(csfit)) {
    # for MCMC samplers
    n_chains <- max(
      nchains(csfit$warmup_draws),
      nchains(csfit$post_warmup_draws)
    )
    n_iter_warmup <- niterations(csfit$warmup_draws)
    n_iter_sample <- niterations(csfit$post_warmup_draws)
    if (n_iter_warmup > 0) {
      csfit$warmup_draws <- as_draws_df(csfit$warmup_draws)
      csfit$warmup_sampler_diagnostics <-
        as_draws_df(csfit$warmup_sampler_diagnostics)
    }
    if (n_iter_sample > 0) {
      csfit$post_warmup_draws <- as_draws_df(csfit$post_warmup_draws)
      csfit$post_warmup_sampler_diagnostics <-
        as_draws_df(csfit$post_warmup_sampler_diagnostics)
    }
    
    # called 'samples' for consistency with rstan
    samples <- rbind(csfit$warmup_draws, csfit$post_warmup_draws)
    # manage memory
    csfit$warmup_draws <- NULL
    csfit$post_warmup_draws <- NULL
    
    # prepare sampler diagnostics
    diagnostics <- rbind(csfit$warmup_sampler_diagnostics,
                         csfit$post_warmup_sampler_diagnostics)
    # manage memory
    csfit$warmup_sampler_diagnostics <- NULL
    csfit$post_warmup_sampler_diagnostics <- NULL
    # convert to regular data.frame
    diagnostics <- as.data.frame(diagnostics)
    diag_chain_ids <- diagnostics$.chain
    diagnostics[res_vars] <- NULL
    
  } else if ("draws" %in% names(csfit)) {
    # for variational inference "samplers"
    n_chains <- 1
    n_iter_warmup <- 0
    n_iter_sample <- niterations(csfit$draws)
    if (n_iter_sample > 0) {
      csfit$draws <- as_draws_df(csfit$draws)
    }
    
    # called 'samples' for consistency with rstan
    samples <- csfit$draws
    # manage memory
    csfit$draws <- NULL
    
    # VI has no sampler diagnostics
    diag_chain_ids <- rep(1L, nrow(samples))
    diagnostics <- as.data.frame(matrix(nrow = nrow(samples), ncol = 0))
  }
  
  # convert to regular data.frame
  samples <- as.data.frame(samples)
  chain_ids <- samples$.chain
  samples[res_vars] <- NULL
  if ("lp__" %in% colnames(samples)) {
    samples <- move2end(samples, "lp__")
  }
  
  fnames_oi <- colnames(samples)
  
  colnames(samples) <- gsub("\\[", ".", colnames(samples))
  colnames(samples) <- gsub("\\]", "", colnames(samples))
  colnames(samples) <- gsub("\\,", ".", colnames(samples))
  
  # split samples into chains
  samples <- split(samples, chain_ids)
  names(samples) <- NULL
  
  # split diagnostics into chains
  diagnostics <- split(diagnostics, diag_chain_ids)
  names(diagnostics) <- NULL
  
  #  @sim$sample: largely 113-130 from rstan::read_stan_csv
  values <- list()
  values$algorithm <- csfit$metadata$algorithm
  values$engine <- csfit$metadata$engine
  values$metric <- csfit$metadata$metric
  
  sampler_t <- NULL
  if (!is.null(values$algorithm)) {
    if (values$algorithm == "rwm" || values$algorithm == "Metropolis") {
      sampler_t <- "Metropolis"
    } else if (values$algorithm == "hmc") {
      if (values$engine == "static") {
        sampler_t <- "HMC"
      } else {
        if (values$metric == "unit_e") {
          sampler_t <- "NUTS(unit_e)"
        } else if (values$metric == "diag_e") {
          sampler_t <- "NUTS(diag_e)"
        } else if (values$metric == "dense_e") {
          sampler_t <- "NUTS(dense_e)"
        }
      }
    }
  }
  
  adapt_info <- vector("list", 4)
  idx_samples <- (n_iter_warmup + 1):(n_iter_warmup + n_iter_sample)
  
  for (i in seq_along(samples)) {
    m <- colMeans(samples[[i]][idx_samples, , drop=FALSE])
    rownames(samples[[i]]) <- seq_rows(samples[[i]])
    attr(samples[[i]], "sampler_params") <- diagnostics[[i]][rstan_diagn_order]
    rownames(attr(samples[[i]], "sampler_params")) <- seq_rows(diagnostics[[i]])
    
    # reformat back to text
    if (is_equal(sampler_t, "NUTS(dense_e)")) {
      mmatrix_txt <- "\n# Elements of inverse mass matrix:\n# "
      mmat <- paste0(apply(csfit$inv_metric[[i]], 1, paste0, collapse=", "),
                     collapse="\n# ")
    } else {
      mmatrix_txt <- "\n# Diagonal elements of inverse mass matrix:\n# "
      mmat <- paste0(csfit$inv_metric[[i]], collapse = ", ")
    }
    
    adapt_info[[i]] <- paste0("# Step size = ",
                              csfit$step_size[[i]],
                              mmatrix_txt,
                              mmat, "\n# ")
    
    attr(samples[[i]], "adaptation_info") <- adapt_info[[i]]
    
    attr(samples[[i]], "args") <- list(sampler_t = sampler_t, chain_id = i)
    
    if (NROW(csfit$metadata$time)) {
      time_i <- as.double(csfit$metadata$time[i, c("warmup", "sampling")])
      names(time_i) <- c("warmup", "sample")
      attr(samples[[i]], "elapsed_time") <- time_i
    }
    
    attr(samples[[i]], "mean_pars") <- m[-length(m)]
    attr(samples[[i]], "mean_lp__") <- m["lp__"]
  }
  
  perm_lst <- lapply(seq_len(n_chains), function(id) sample.int(n_iter_sample))
  
  # @sim
  sim <- list(
    samples = samples,
    iter = csfit$metadata$iter_sampling + csfit$metadata$iter_warmup,
    thin = csfit$metadata$thin,
    warmup = csfit$metadata$iter_warmup,
    chains = n_chains,
    n_save = rep(n_iter_sample + n_iter_warmup, n_chains),
    warmup2 = rep(n_iter_warmup, n_chains),
    permutation = perm_lst,
    pars_oi = pars_oi,
    dims_oi = par_dims,
    fnames_oi = fnames_oi,
    n_flatnames = length(fnames_oi)
  )
  
  # @stan_args
  sargs <- list(
    stan_version_major = as.character(csfit$metadata$stan_version_major),
    stan_version_minor = as.character(csfit$metadata$stan_version_minor),
    stan_version_patch = as.character(csfit$metadata$stan_version_patch),
    model = csfit$metadata$model_name,
    start_datetime = gsub(" ", "", csfit$metadata$start_datetime),
    method = csfit$metadata$method,
    iter = csfit$metadata$iter_sampling + csfit$metadata$iter_warmup,
    warmup = csfit$metadata$iter_warmup,
    save_warmup = csfit$metadata$save_warmup,
    thin = csfit$metadata$thin,
    engaged = as.character(csfit$metadata$adapt_engaged),
    gamma = csfit$metadata$gamma,
    delta = csfit$metadata$adapt_delta,
    kappa = csfit$metadata$kappa,
    t0 = csfit$metadata$t0,
    init_buffer = as.character(csfit$metadata$init_buffer),
    term_buffer = as.character(csfit$metadata$term_buffer),
    window = as.character(csfit$metadata$window),
    algorithm = csfit$metadata$algorithm,
    engine = csfit$metadata$engine,
    max_depth = csfit$metadata$max_treedepth,
    metric = csfit$metadata$metric,
    metric_file = character(0), # not stored in metadata
    stepsize = NA, # add in loop
    stepsize_jitter = csfit$metadata$stepsize_jitter,
    num_chains = as.character(csfit$metadata$num_chains),
    chain_id = NA, # add in loop
    file = character(0), # not stored in metadata
    init = NA, # add in loop
    seed = as.character(csfit$metadata$seed),
    file = NA, # add in loop
    diagnostic_file = character(0), # not stored in metadata
    refresh = as.character(csfit$metadata$refresh),
    sig_figs = as.character(csfit$metadata$sig_figs),
    profile_file = csfit$metadata$profile_file,
    num_threads = as.character(csfit$metadata$threads_per_chain),
    stanc_version = gsub(" ", "", csfit$metadata$stanc_version),
    stancflags = character(0), # not stored in metadata
    adaptation_info = NA, # add in loop
    has_time = is.numeric(csfit$metadata$time$total),
    time_info = NA, # add in loop
    sampler_t = sampler_t
  )
  
  sargs_rep <- replicate(n_chains, sargs, simplify = FALSE)
  
  for (i in seq_along(sargs_rep)) {
    sargs_rep[[i]]$chain_id <- i
    sargs_rep[[i]]$stepsize <- csfit$metadata$step_size[i]
    sargs_rep[[i]]$init <- as.character(csfit$metadata$init[i])
    # two 'file' elements: select the second
    file_idx <- which(names(sargs_rep[[i]]) == "file")
    sargs_rep[[i]][[file_idx[2]]] <- files[[i]]
    
    sargs_rep[[i]]$adaptation_info <- adapt_info[[i]]
    
    if (NROW(csfit$metadata$time)) {
      sargs_rep[[i]]$time_info <- paste0(
        c("#  Elapsed Time: ", "#                ", "#                ", "# "),
        c(csfit$metadata$time[i, c("warmup", "sampling", "total")], ""),
        c(" seconds (Warm-up)", " seconds (Sampling)", " seconds (Total)", "")
      )
    }
  }
  
  # @stanmodel
  null_dso <- new(
    "cxxdso", sig = list(character(0)), dso_saved = FALSE,
    dso_filename = character(0), modulename = character(0),
    system = R.version$system, cxxflags = character(0),
    .CXXDSOMISC = new.env(parent = emptyenv())
  )
  null_sm <- new(
    "stanmodel", model_name = model_name, model_code = character(0),
    model_cpp = list(), dso = null_dso
  )
  
  # @date
  sdate <- do.call(max, lapply(files, function(csv) file.info(csv)$mtime))
  sdate <- format(sdate, "%a %b %d %X %Y")
  
  new(
    "stanfit",
    model_name = model_name,
    model_pars = svars,
    par_dims = par_dims,
    mode = mode,
    sim = sim,
    inits = list(),
    stan_args = sargs_rep,
    stanmodel = null_sm,
    date = sdate,  # not the time of sampling
    .MISC = new.env(parent = emptyenv())
  )
}

ulapply <- function(X, FUN, ..., recursive = TRUE, use.names = TRUE) {
  unlist(lapply(X, FUN, ...), recursive, use.names)
}







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
