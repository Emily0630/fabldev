
multivar_beta <- function(x){
    sum(lgamma(x)) - lgamma(sum(x))
}
samps <- 50
alpha <- rep(1, ncol(hash$ohe))
beta <- rep(1, ncol(hash$ohe))
alpha_pi <- 1
beta_pi <- 1
phis <- c(out$pattern_weights, exp(digamma(out$b_pi)))
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

set.seed(42)
log_ratios <- rep(0, 10000)

for(i in 1:10000){
  # fabl parameters
  m <- as.vector(unlist(sapply(a_split, function(x){
    prob <- MCMCpack::rdirichlet(1, x)
    prob/sum(prob)
  })))

  u <- as.vector(unlist(sapply(b_split, function(x){
    prob <- MCMCpack::rdirichlet(1, x)
    prob/sum(prob)
  })))

  pi <- rbeta(1, a_pi, b_pi)

  log_m_p <- ohe %>%
    sweep(., 2, log(m), "*") %>%
    rowSums()

  log_u_p <- ohe %>%
    sweep(., 2, log(u), "*") %>%
    rowSums()


  xi <- sapply(1:n2, function(j){
    sample(1:(P+1), 1, F, pattern_probs[[j]])
  })

  n_p_xi <- xi %>%
    factor(., 1:(P+1), 1:(P+1)) %>%
    table() %>%
    .[-(P+1)] %>%
    as.numeric()

  n_12_xi <- sum(n_p_xi)

  #importance ratio
  log_ratio <- sum(n_p_xi * log_m_p) +
    sum((N_p - n_p_xi) * log_u_p) +
    sum((alpha - a) * log(m)) +
    sum((beta - b) * log(u)) +
    sum(sapply(a_split, multivar_beta)) +
    sum(sapply(b_split, multivar_beta)) -
    sum(sapply(alpha_split, multivar_beta)) -
    sum(sapply(beta_split, multivar_beta)) +
    (alpha_pi - a_pi + n_12_xi) * log(pi) +
    (beta_pi - b_pi + n2 - n_12_xi)  * log(1 - pi) +
    lbeta(a_pi, b_pi) -
    lbeta(alpha_pi, beta_pi) -
    n_12_xi * log(n1) -
    sum(log(phis)[xi] - log(out$C)) #x_j = 0 is coded as P+1

  ratio <- exp(log_ratio)

  log_joint <- sum(n_p_xi * log_m_p) +
    sum((N_p - n_p_xi) * log_u_p) +
    sum((alpha - 1) * log(m)) +
    sum((beta - 1) * log(u)) -
    sum(sapply(alpha_split, multivar_beta)) -
    sum(sapply(beta_split, multivar_beta)) +
    (alpha_pi - 1 + n_12_xi) * log(pi) +
    (beta_pi - 1 + n2 - n_12_xi)  * log(1 - pi) -
    lbeta(alpha_pi, beta_pi) -
    n_12_xi * log(n1)

  log_var <- sum((a - 1) * log(m)) +
    sum((b - 1) * log(u)) -
    sum(sapply(a_split, multivar_beta)) -
    sum(sapply(b_split, multivar_beta)) +
    (a_pi - 1) * log(pi) +
    (b_pi - 1) * log(1 - pi) -
    lbeta(a_pi, b_pi) +
    (n2 - n_12_xi) * digamma(b_pi) +
    sum(n_p_xi * log(out$pattern_weights)) -
    sum(log(out$C))

  log_ratios[i] <- log_joint - log_var
}

out$elbo_seq
mean(log_ratios)


# m <- a %>%
#   split(., field_marker) %>%
#   lapply(., function(x){
#     x/sum(x)
#   }) %>%
#   unlist() %>%
#   log()
#
# u <- b %>%
#   split(., field_marker)%>%
#   lapply(., function(x){
#     x/sum(x)
#   }) %>%
#   unlist() %>%
#   log()

# xi_samps <- lapply(seq_len(samps), function(x){
#   xi <- sapply(1:n2, function(j){
#     sample(1:(P+1), 1, F, pattern_probs[[j]])
#   })
#   xi
# })
#
#
# sapply(list(a, b), function(y){
#   split(y, field_marker) %>%
#     sapply(., function(x){
#       sum(lgamma(x)) - lgamma(sum(x))
#     })%>%
#     sum(.)
# }) %>%
#   sum(.) -
#
#   sapply(list(alpha, beta), function(y){
#     split(y, field_marker) %>%
#       sapply(., function(x){
#         sum(lgamma(x)) - lgamma(sum(x))
#       })%>%
#       sum(.)
#   }) %>%
#   sum(.) +
