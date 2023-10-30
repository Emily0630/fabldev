
samps <- 50
alpha <- rep(1, ncol(hash$ohe))
beta <- rep(1, ncol(hash$ohe))
alpha_pi <- 1
beta_pi <- 1
phis <- c(out$pattern_weights, exp(digamma(out$b_pi)))
n2 <- hash$n2
n1 <- hash$n1
P <- nrow(hash$ohe)
pattern_probs <- lapply(1:n2, function(j){
  c(out$pattern_weights * hash$pattern_counts_by_record[[j]],
    exp(digamma(out$b_pi))) / out$C[j]
})



N_p <- hash$total_counts

a <- out$a
b <- out$b
a_pi <- out$a_pi
b_pi <- out$b_pi
pi <- a_pi / hash$n2
field_marker <- hash$field_marker
ohe <- hash$ohe

# fabl parameters

m <- a %>%
  split(., field_marker) %>%
  lapply(., function(x){
    x/sum(x)
  }) %>%
  unlist() %>%
  log()

u <- b %>%
  split(., field_marker)%>%
  lapply(., function(x){
    x/sum(x)
  }) %>%
  unlist() %>%
  log()

m_p <- ohe %>%
  sweep(., 2, m, "*") %>%
  rowSums() %>%
  exp()


u_p <- ohe %>%
  sweep(., 2, u, "*") %>%
  rowSums() %>%
  exp()

# xi

xi_samps <- lapply(seq_len(samps), function(x){
  xi <- sapply(1:n2, function(j){
    sample(1:(P+1), 1, F, pattern_probs[[j]])
  })
  xi
})

xi <- sapply(1:n2, function(j){
  sample(1:(P+1), 1, F, pattern_probs[[j]])
})
xi

n_p_xi <- xi %>%
  factor(., 1:(P+1), 1:(P+1)) %>%
  table() %>%
  .[-(P+1)] %>%
  as.numeric()
n_12_xi <- sum(n_p_xi)

#importance ratio

log_ratio <- sum(n_p_xi * log(m_p) + (N_p - n_p_xi) * log(u_p)) +
  sum((alpha - a) * m + (beta - b) * u) +
  sapply(list(a, b), function(y){
  split(y, field_marker) %>%
    sapply(., function(x){
      sum(lgamma(x)) - lgamma(sum(x))
    })%>%
    sum(.)
}) %>%
  sum(.) -
  sapply(list(alpha, beta), function(y){
  split(y, field_marker) %>%
    sapply(., function(x){
      sum(lgamma(x)) - lgamma(sum(x))
    })%>%
    sum(.)
}) %>%
  sum(.) +
  (alpha_pi - a_pi + n_12_xi) * log(pi) +
  (beta_pi - b_pi + n2 - n_12_xi)  * log(1 - pi) +
  lbeta(a_pi, b_pi) -
  lbeta(alpha_pi, beta_pi) -
  n_12_xi * log(n1) - sum(log(phis)[xi] / out$C)

ratio <- exp(log_ratio)
