stan_hurdle_cumulative_lpmf <- function(family, link) {
  stopifnot(is.character(family), is.character(link))
  inv_link <- stan_inv_link(link)
  th <- function(k) {
    out <- glue("thres[{k}] - mu")
    glue("disc * ({out})")
  }
  out <- glue(
    "  /* {family}-{link} log-PDF for a single response\n",
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: latent mean parameter\n",
    "   *   hu: hurdle probability\n",
    "   *   disc: discrimination parameter\n",
    "   *   thres: ordinal thresholds\n",
    "   * Returns:\n",
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}_lpmf(int y, real mu, real hu, real disc, vector thres) {{\n",
    "\n"
  ) 
  # define the function body
  if (inv_link == "inv_logit") {
    str_add(out) <- glue(
      "     int nthres = num_elements(thres);\n",
      "     if (y == 0) {{\n",
      "       return bernoulli(1 | hu);\n",
      "     }} else if (y == 1) {{\n",
      "       return log_inv_logit({th(1)}) +\n", 
      "                bernoulli(0 | hu);\n",
      "     }} else if (y == nthres + 2) {{\n",
      "       return log1m_inv_logit({th('nthres')}) +\n",
      "                bernoulli(0 | hu);\n",
      "     }} else {{\n",
      # TODO: replace with log_inv_logit_diff once rstan >= 2.25
      "       return log_diff_exp(\n",
      "         log_inv_logit({th('y')}), \n",
      "         log_inv_logit({th('y - 1')})\n",
      "       ) + bernoulli(0 | hu) ;\n",
      "     }}\n",
      "   }}\n"
    )
  } else {
    str_add(out) <- glue(
      "     int nthres = num_elements(thres);\n",
      "     real p;\n",
      "     if (y == 0){{\n",
      "       p = hu;\n",
      "     else if (y == 1) {{\n",
      "       p = {inv_link}({th(1)}) * (1 - hu);\n",
      "     }} else if (y == nthres + 2) {{\n",
      "       p = (1 - {inv_link}({th('nthres')})) * (1 - hu);\n",
      "     }} else {{\n",
      "       p = ({inv_link}({th('y')}) -\n",
      "           {inv_link}({th('y - 1')})) * (1 - hu);\n",
      "     }}\n",
      "     return log(p);\n",
      "   }}\n"
    )
    
  } 
  
  # Use more efficient ordered_logistic function when disc == 1
  str_add(out) <- glue(
        "   real hurdle_cumulative_ordered_logistic_lpmf(int y, real mu, real hu, real disc, vector thres) {{\n",
        "     if (y == 0) {{\n",
        "       return bernoulli(1 | hu);\n",
        "     }} else {{\n",
        "       return ordered_logistic_lpmf(y | mu, thres) +\n", 
        "                bernoulli(0 | hu);\n",
        "     }}\n",
        "   }}\n"
      )
  
  # lpmf function for multiple merged thresholds
  str_add(out) <- glue(
    "  /* {family}-{link} log-PDF for a single response and merged thresholds\n",
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: latent mean parameter\n",
    "   *   hu: hurdle probability\n",
    "   *   disc: discrimination parameter\n",
    "   *   thres: vector of merged ordinal thresholds\n",
    "   *   j: start and end index for the applid threshold within 'thres'\n",
    "   * Returns:\n",
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}_merged_lpmf(",
    "int y, real mu, real hu, real disc, vector thres, int[] j) {{\n",
    "     return {family}_{link}_lpmf(y | mu, hu, disc, thres[j[1]:j[2]]);\n",
    "   }}\n"
  )
  
  if (link == "logit") {
    # use the more efficient 'ordered_logistic' built-in function
    str_add(out) <- glue(
      "  /* use ordered-logistic log-PDF for a single response and merged thresholds\n",
      "   * Args:\n",
      "   *   y: response category\n",
      "   *   mu: latent mean parameter\n",
      "   *   hu: hurdle probability\n",
      "   *   thres: vector of merged ordinal thresholds\n",
      "   *   j: start and end index for the applid threshold within 'thres'\n",
      "   * Returns:\n",
      "   *   a scalar to be added to the log posterior\n",
      "   */\n",
      "   real hurdle_cumulative_ordered_logistic_merged_lpmf(",
      "int y, real mu, real hu, vector thres, int[] j) {{\n",
      "     return hurdle_cumulative_ordered_logistic_lpmf(y | mu, hu, thres[j[1]:j[2]]);\n",
      "   }}\n"
    )
  }
  out
}


else {
  str_add(out) <- glue(
    "     int nthres = num_elements(thres);\n",
    "     real p;\n",
    "     if (y == 1) {{\n",
    "       p = {inv_link}({th(1)});\n",
    "     }} else if (y == nthres + 1) {{\n",
    "       p = 1 - {inv_link}({th('nthres')});\n",
    "     }} else {{\n",
    "       p = {inv_link}({th('y')}) -\n",
    "           {inv_link}({th('y - 1')});\n",
    "     }}\n",
    "     return log(p);\n",
    "   }}\n"
  )
}

# lpmf function for multiple merged thresholds
str_add(out) <- glue(
  "  /* {family}-{link} log-PDF for a single response and merged thresholds\n",
  "   * Args:\n",
  "   *   y: response category\n",
  "   *   mu: latent mean parameter\n",
  "   *   disc: discrimination parameter\n",
  "   *   thres: vector of merged ordinal thresholds\n",
  "   *   j: start and end index for the applid threshold within 'thres'\n",
  "   * Returns:\n",
  "   *   a scalar to be added to the log posterior\n",
  "   */\n",
  "   real {family}_{link}_merged_lpmf(",
  "int y, real mu, real disc, vector thres, int[] j) {{\n",
  "     return {family}_{link}_lpmf(y | mu, disc, thres[j[1]:j[2]]);\n",
  "   }}\n"
)
if (family == "cumulative" && link == "logit") {
  # use the more efficient 'ordered_logistic' built-in function
  str_add(out) <- glue(
    "  /* ordered-logistic log-PDF for a single response and merged thresholds\n",
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: latent mean parameter\n",
    "   *   thres: vector of merged ordinal thresholds\n",
    "   *   j: start and end index for the applid threshold within 'thres'\n",
    "   * Returns:\n",
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real ordered_logistic_merged_lpmf(",
    "int y, real mu, vector thres, int[] j) {{\n",
    "     return ordered_logistic_lpmf(y | mu, thres[j[1]:j[2]]);\n",
    "   }}\n"
  )
}
out
}




d <- data.frame(y = c(0, 1),
                x = c(0, 1),
                ntrials = c(5, 5))

make_stancode(y ~ x, family = bernoulli(link = "cauchit"), data = d)
make_stancode(bf(y | Trials(ntrials) ~ x, zi ~ x), family = zero_inflated_binomial(link_zi = "cauchit"), data = d)







stan_inv_link <- function(link) {
  switch(link,
         identity = "",
         log = "exp",
         logm1 = "expp1",
         inverse = "inv",
         sqrt = "square",
         "1/mu^2" = "inv_sqrt",
         logit = "inv_logit",
         probit = "Phi",
         probit_approx = "Phi_approx",
         cloglog = "inv_cloglog",
         cauchit = "inv_cauchit",
         tan_half = "inv_tan_half",
         log1p = "expm1",
         softplus = "log1p_exp",
         squareplus = "squareplus",
         softit = "inv_softit"
  )
}

'str_add<-' <- function(x, start = FALSE, value) {
  if (start) paste0(value, x) else paste0(x, value)
}
