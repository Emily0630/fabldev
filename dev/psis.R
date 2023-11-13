multivar_beta <- function(x){
    sum(lgamma(x)) - lgamma(sum(x))
}

psis_vabl <- function(hash, out, num_samp, seed = 42){

  # Get prior information and variational estimates
  alpha <- rep(1, ncol(hash$ohe))
  beta <- rep(1, ncol(hash$ohe))
  alpha_pi <- 1
  beta_pi <- 1
  phis <- c(out$pattern_weights, exp(digamma(out$b_pi)))
  log_phis <- log(phis)
  n2 <- hash$n2
  n1 <- hash$n1
  P <- nrow(hash$ohe)
  field_marker <- hash$field_marker
  ohe <- hash$ohe
  pattern_probs <- lapply(1:n2, function(j){
    c(out$pattern_weights * hash$pattern_counts_by_record[[j]],
      exp(digamma(out$b_pi))) / out$C[j]
  })

  N_p <- hash$total_counts

  a <- out$a
  b <- out$b
  a_split <- split(a, field_marker)
  b_split <- split(b, field_marker)
  a_pi <- out$a_pi
  b_pi <- out$b_pi
  alpha_split <- split(alpha, field_marker)
  beta_split <- split(beta, field_marker)

  # Pre-compute reusable things
  beta_funs_mu <- sum(sapply(a_split, multivar_beta)) +
    sum(sapply(b_split, multivar_beta)) -
    sum(sapply(alpha_split, multivar_beta)) -
    sum(sapply(beta_split, multivar_beta))

  beta_funs_pi <- lbeta(a_pi, b_pi) -
    lbeta(alpha_pi, beta_pi)

  Phi_sum <- sum(log(out$C))

  set.seed(seed)

  log_ms <- log(replicate(num_samp, as.vector(unlist(sapply(a_split, function(x){
    MCMCpack::rdirichlet(1, x)
  })))))
  log_us <- log(replicate(num_samp, as.vector(unlist(sapply(b_split, function(x){
    MCMCpack::rdirichlet(1, x)
  })))))
  log_m_ps <- ohe %*% log_ms
  log_u_ps <- ohe %*% log_us

  pis <- rbeta(num_samp, a_pi, b_pi)

  xis <- replicate(num_samp, sapply(1:n2, function(j){
    sample(1:(P+1), 1, F, pattern_probs[[j]])
  }))
  n_p_xis <- apply(xis, 2, function(x){
    tabulate(bin = x, nbin = P + 1)[-(P + 1)]
  })

  log_ratios <- rep(0, num_samp)

  for(i in 1:num_samp){
    # Sample parameters

    # m <- as.vector(unlist(sapply(a_split, function(x){
    #   MCMCpack::rdirichlet(1, x)
    # })))
    #
    # u <- as.vector(unlist(sapply(b_split, function(x){
    #   MCMCpack::rdirichlet(1, x)
    # })))

    log_m <- log_ms[, i]
    log_u <- log_us[, i]

    log_m_p <- log_m_ps[, i]
    log_u_p <- log_u_ps[, i]

    pi <- pis[i]

    # xi <- sapply(1:n2, function(j){
    #   sample(1:(P+1), 1, F, pattern_probs[[j]])
    # })

    xi <- xis[, i]
    # n_p_xi <- xi %>%
    #   factor(., 1:(P+1), 1:(P+1)) %>%
    #   table() %>%
    #   .[-(P+1)] %>%
    #   as.numeric()
    n_p_xi <- n_p_xis[, i]

    n_12_xi <- sum(n_p_xi)

    # Compute log ratio
    # log_ratios[i] <- sum(n_p_xi * log_m_p) +
    #   sum((N_p - n_p_xi) * log_u_p) +
    #   sum((alpha - a) * log(m)) +
    #   sum((beta - b) * log(u)) +
    #   # beta_funs_mu +
    #   sum(sapply(a_split, multivar_beta)) +
    #   sum(sapply(b_split, multivar_beta)) -
    #   sum(sapply(alpha_split, multivar_beta)) -
    #   sum(sapply(beta_split, multivar_beta)) +
    #   (alpha_pi - a_pi + n_12_xi) * log(pi) +
    #   (beta_pi - b_pi + n2 - n_12_xi)  * log(1 - pi) +
    #   lbeta(a_pi, b_pi) -
    #   lbeta(alpha_pi, beta_pi) -
    #   # beta_funs_pi -
    #   n_12_xi * log(n1) -
    #   sum(log(phis)[xi]) -
    #   sum(log(out$C))
    #   # Phi_sum
    log_ratios[i] <- sum(n_p_xi * log_m_p) +
      sum((N_p - n_p_xi) * log_u_p) +
      sum((alpha - a) * log_m) +
      sum((beta - b) * log_u) +
      beta_funs_mu +
      (alpha_pi - a_pi + n_12_xi) * log(pi) +
      (beta_pi - b_pi + n2 - n_12_xi)  * log(1 - pi) +
      beta_funs_pi -
      n_12_xi * log(n1) -
      sum(log_phis[xi]) + Phi_sum
  }

  return(log_ratios)
}


num_samp <- 10000
log_ratios <- psis_vabl(hash = hash, out = out, num_samp = num_samp, seed = 43)
hist(log_ratios)
hist(log_ratios - max(log_ratios))
out$elbo_seq
mean(log_ratios)
# diagnostic <- loo::psis(log_ratios = log_ratios, r_eff = NA)
# diagnostic
# diagnostic$diagnostic$pareto_k

# Was getting weird results for febrl when using the psis function
# I know the issue (has to do with a cutoff), so using an older version that
# uses a different cutoff
M <- ceiling(min(num_samp / 5, 3 * sqrt(num_samp)))
k <- loo:::psislw(log_ratios, wcp = M / num_samp)$pareto_k


log_ratios <- psis_vabl(hash = hash, out = out, num_samp = 100000 * 100, seed = 42)
khats <- rep(0, 20)
for(j in 1:20){
  khats[j] <- loo::psis(log_ratios = log_ratios[(500000 * (j-1) + 1):(500000 * j)],
                        r_eff = NA)$diagnostic$pareto_k
}
hist(khats)
