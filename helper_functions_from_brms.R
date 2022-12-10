

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

logit <- function(x) log(x / (1 - x))

slice_col <- function(x, i) {
  if (length(dim(x)) < 2L) {
    # a vector or 1D array is treated as already sliced
    return(x)
  }
  slice(x, 2, i)
}


slice <- function(x, dim, i, drop = TRUE) {
  ndim <- length(dim(x))
  commas1 <- collapse(rep(", ", dim - 1))
  commas2 <- collapse(rep(", ", ndim - dim))
  drop_dim <- ifelse(drop, ", drop_dim = dim", "")
  expr <- paste0("extract(x, ", commas1, "i", commas2, drop_dim, ")")
  eval2(expr)
}


do_call <- function(what, args, pkg = NULL, envir = parent.frame()) {
  call <- ""
  if (length(args)) {
    if (!is.list(args)) {
      stop2("'args' must be a list.")
    }
    fun_args <- names(args)
    if (is.null(fun_args)) {
      fun_args <- rep("", length(args))
    } else {
      nzc <- nzchar(fun_args)
      fun_args[nzc] <- paste0("`", fun_args[nzc], "` = ")
    }
    names(args) <- paste0(".x", seq_along(args))
    call <- paste0(fun_args, names(args), collapse = ",")
  } else {
    args <- list()
  }
  if (is.function(what)) {
    args$.fun <- what
    what <- ".fun"
  } else {
    what <- paste0("`", as_one_character(what), "`")
    if (!is.null(pkg)) {
      what <- paste0(as_one_character(pkg), "::", what)
    }
  }
  call <- paste0(what, "(", call, ")")
  eval2(call, envir = args, enclos = envir)
}

eval2 <- function(expr, envir = parent.frame(), ...) {
  if (is.character(expr)) {
    expr <- str2expression(expr)
  }
  eval(expr, envir, ...)
}
