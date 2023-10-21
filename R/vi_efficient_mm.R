vi_efficient_mm <- function(hash, threshold = 1e-5, tmax = 1000,
                            fixed_iterations = NULL,
                           b_init = TRUE, use_elbo = F){

  check_every <- 10

  ohe <- hash$ohe # One-hot encodings e(h_p)
  combo_ohe <- hash$combination_ohe
  combination_counts <- hash$combination_counts

  P <- dim(ohe)[1]
  P_vec <- sapply(combo_ohe, nrow)
  total_counts <- hash$total_counts #N_p
  pattern_counts_by_record <- hash$pattern_counts_by_record #N_p_j
  record_counts_by_pattern <- hash$record_counts_by_pattern
  field_marker <- hash$field_marker
  pattern_finder <- hash$pattern_finder
  n1 <- hash$n1
  n2 <- hash$n2
  max_K <- hash$max_K
  n1_vec <- sapply(seq_len(max_K), function(x){
    choose(n1, x)
  })


  # Priors
  alpha <- rep(1, length(field_marker))
  Beta <- rep(1, length(field_marker))
  alpha_pi <- rep(1, max_K)
  beta_pi <- rep(1, max_K)
  # alpha_pi <- rep(1, max_K+1)
  # beta_pi <- rep(1, max_K+1)

  # Initialize
    a <- rep(1, length(field_marker))
  if(b_init == T){
    b <- hash$ohe %>%
      sweep(., 1, hash$total_counts, "*") %>%
      colSums()
  } else {
    b = rep(1, length(field_marker))
  }
    a_pi <- rep(1, max_K)
    b_pi <- rep(1, max_K)
    # a_pi <- rep(1, max_K +1)
    # b_pi <- rep(1, max_K+1)
    total_nonmatch <- n2
    total_match <- 0
    #b_pi[max_K] <- NA
    #a_pi <- rep(1, max_K + 1)
    #a_pi <- c(1, rep(1/max_K, max_K))

  t <- 1
  ratio <- 1
  elbo_seq <- vector()

  while(t <= tmax){
    a_sum <- a %>%
      split(., field_marker) %>%
      sapply(., sum) %>%
      digamma(.) %>%
      .[field_marker]

    a_chunk <- digamma(a) - a_sum

    b_sum <- b %>%
      split(., field_marker) %>%
      sapply(., sum) %>%
      digamma(.) %>%
      .[field_marker]
    b_chunk <- digamma(b) - b_sum

    m_p <- lapply(combo_ohe, function(x){
      sweep(x, 2, a_chunk, "*") %>%
      rowSums()
    })

    u_p <- lapply(combo_ohe, function(x){
      sweep(x, 2, b_chunk, "*") %>%
        rowSums()
    })

    # w_p
    weights = purrr::map2(m_p, u_p, ~ .x -.y)

    phi <- sapply(seq_len(max_K), function(k){
      exp(digamma(a_pi[k]) - digamma(n1_vec[k]) + weights[[k]]
          -digamma(b_pi[k]))
    })
    # phi <- sapply(seq_len(max_K), function(k){
    #   exp(digamma(a_pi[k+1]) - digamma(n1_vec[k]) + weights[[k]]
    #       -digamma(b_pi[k+1]))
    # })
    #single <- exp(digamma(total_nonmatch + 1) - digamma(n2 - total_nonmatch + 1))
    #single <- exp(digamma(a_pi[1]) - digamma(b_pi[1]))
    #single <- exp(digamma(total_nonmatch + 1) - digamma(total_match + 1))
    single <- 1

    # Phi_j
    C_temp <- sapply(seq_len(max_K), function(k){
      sapply(combination_counts, function(x){
        x[[k]] %*% phi[[k]]
      })
    })

    C <- C_temp%>%
      rowSums() + single

    # S(Phi)
    temp_match <- colSums(cbind(unname(single), C_temp) / C)
    #cumulative_match <- cumsum(temp_match[max_K:1])[max_K:1]
    cumulative_match <- cumsum(temp_match[(max_K+1):1])[(max_K+1):1]
    a_mass <- cumulative_match[-1]
    b_mass <- cumulative_match[-(max_K+1)] - a_mass
    total_match <- cumulative_match[1]
    total_nonmatch <- n2 - cumulative_match[1]

    # N_p(Psi)
    K_big <- lapply(seq_len(max_K), function(k){
      sapply(seq_len(P_vec[k]), function(p){
        sum(record_counts_by_pattern[[k]][[p]]/C)
      })
    }) %>%
      unlist()

    phi_big <- unlist(phi)

    proper_weight <- K_big * phi_big

    K_proper <- sapply(1:P_vec[1], function(y){
      sum(proper_weight[pattern_finder[[y]]])
    })

    # thing <- combination_counts[[6]] %>% unlist()
    # max(phi_big * thing/ C[6])
    #

    AZ <- ohe %>%
      sweep(., 1, K_proper, "*") %>%
      colSums()

    BZ <- ohe %>%
      sweep(., 1, total_counts - K_proper, "*") %>%
      colSums()

    a <- alpha + AZ
    b <- Beta+ BZ

    #a_pi <- alpha_pi + cumulative_match
    #b_pi <- beta_pi + c(cumulative_match[-1], 0)
    a_pi <- alpha_pi + a_mass
    b_pi <- beta_pi + b_mass
    #b_pi[max_K] <- NA
    t = t + 1

}
    # ELBO

    if (use_elbo == T){
    elbo_pieces <- vector(length = 6)
    elbo_pieces[1] <- sapply(1:n2, function(j){
      sum(pattern_counts_by_record[[j]] * (phi *(weights - log(phi) + log(C[j]))/ C[j] +
                                             u_p))
    }) %>%
      sum(.)

    #elbo_pieces[2] <- -sum(single/C *log(single/C)) + total_nonmatch * log(n1) -log(n1)*n2
    elbo_pieces[2] <- single * sum(1/C *log(C)) + total_nonmatch * (log(n1) - log(single)) -log(n1)*n2

    elbo_pieces[3] <- lbeta(a_pi, b_pi) - lbeta(alpha_pi, beta_pi)


    elbo_pieces[4] <- sapply(list(a, b), function(y){
      split(y, field_marker) %>%
        sapply(., function(x){
          sum(lgamma(x)) - lgamma(sum(x))
        })%>%
        sum(.)
    }) %>%
      sum(.)

    elbo_pieces[5] <- -sapply(list(alpha, Beta), function(y){
      split(y, field_marker) %>%
        sapply(., function(x){
          sum(lgamma(x)) - lgamma(sum(x))
        })%>%
        sum(.)
    }) %>%
      sum(.)

    elbo_pieces[6] <- sum((alpha - a) * a_chunk + (Beta - b) * b_chunk)
    elbo <- sum(elbo_pieces)
    elbo_seq <- c(elbo_seq, elbo)


    if(is.null(fixed_iterations)){
      if(t %% check_every == 0){
        ratio <- abs((elbo_seq[t] - elbo_seq[t - check_every +1])/
                       elbo_seq[t - check_every +1])
      }
      if(ratio < threshold){
        break
      }
    }

    t <- t + 1
    if(t > tmax){
      print("Max iterations have passed before convergence")
      break
    }

    if(!is.null(fixed_iterations)){
      if(t == fixed_iterations){
        break
      }
    }


  }

  list(pattern_weights = phi_big,
       C = C,
       a = a,
       b = b,
       a_pi = a_pi,
       b_pi = b_pi,
       elbo_seq = elbo_seq,
       t = t)

}

# record <- 11
# thing <- combination_counts[[record]] %>% unlist()
# probs <- phi_big * thing/ C[record]
# probs_reduced <- sapply(1:P, function(x){
#   sum(probs[pattern_finder[[x]]])
# })
# match_class <- which(probs_reduced > .5)
#
# matching_records <- lapply(match_class, function(x){
#   hash$hash_to_file_1[[record]][[x]]
# }) %>%
#   unlist()
# matching_records
