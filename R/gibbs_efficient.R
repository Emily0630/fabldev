#' @export
#'

gibbs_efficient <- function(hash, m_prior = 1, u_prior = 1,
                                  alpha = 1, beta = 1, S = 1000, burn = round(S * .1),
                                  show_progress = T){
  # Implements bipartite record linkage with BK Sampling Mechanism
  #
  # Arguments
  # comparisons = list calculated from from BRL::compareRecords
  # m.prior = prior distribution for m parameters
  # u.prior= prior distribution for u parameters
  # alpha = first parameter of prior for linkage probability
  # beta = second parameter of prior for linkage probability
  # S = number of Gibbs iterations
  # burn = number of iterations to be discarded as burn-in
  # show_progress = set to false to show simulation progress

  n1 <- hash$n1
  n2 <- hash$n2
  field_marker <- hash$field_marker

  unique_patterns <- hash$ohe
  pattern_counts <- hash$total_counts
  P <- nrow(unique_patterns)
  hash_count_list <- hash$hash_count_list
  #counts_by_rec <- hash$pattern_counts_by_record
  hash_to_file_1 <-hash$hash_to_file_1

  #candidates_P <- 1:(P+1)
  candidates_P <- 0:P
  Z.SAMPS <- matrix(0, nrow = n2, ncol = S)
  M.SAMPS <- matrix(NA, nrow = length(field_marker), ncol = S)
  U.SAMPS <- matrix(NA, nrow = length(field_marker), ncol = S)
  L.SAMPS <- vector(length = S)
  PI.SAMPS <- vector(length = S)
  Z <- rep(0, n2)
  L <- 0

  m <- u <- rep(0, length(field_marker))
  matches <- rep(0,P)
  #set.seed(1)

  # Gibbs
  for(s in 1:S){

    AZ <- sweep(unique_patterns, MARGIN = 1, STAT = matches, FUN = "*") %>%
      colSums() %>%
      unname()

    nonmatches <- pattern_counts - matches

    BZ <- sweep(unique_patterns, MARGIN = 1, STAT = nonmatches, FUN = "*") %>%
      colSums() %>%
      unname()


    m_post <- m_prior + AZ
    u_post <- u_prior + BZ

    m_post <- split(m_post, field_marker)
    m <- as.vector(unlist(sapply(m_post, function(x){
      prob <- MCMCpack::rdirichlet(1, x)
      prob/sum(prob)
    })))

    u_post <- split(u_post, field_marker)
    u <- as.vector(unlist(sapply(u_post, function(x){
      prob <- MCMCpack::rdirichlet(1, x)
      prob/sum(prob)
    })))


    ratio <- (log(m) - log(u)) %>%
      rep(., P) %>%
      matrix(., nrow = P, byrow = TRUE)

    unique_weights <- exp(rowSums(ratio * unique_patterns, na.rm = TRUE))

    # hash_weights <- apply(hash_count_table, 2, function(x){
    #   x * unique_weights
    # })

    hash_weights <- lapply(hash_count_list, function(x){
      x * unique_weights
    })

    pi <- rbeta(1, L + alpha, n2 - L + beta)

    for(j in 1:n2){
      if(Z[j] > 0){
        L <- L - 1
      }
      Z[j] <- sample(candidates_P, 1,
                     prob = c(1 - pi, hash_weights[[j]] * pi / n1))
      if(Z[j] > 0){
        #index <- ceiling(runif(1) * counts_by_rec[[j]][Z[j]]) #Causes issue with SEI
        index <- ceiling(runif(1) * length(hash_to_file_1[[j]][[Z[j]]]))
        Z.SAMPS[j, s] <- hash_to_file_1[[j]][[Z[j]]][index]
        L <- L + 1
      }
    }
    hash_matches <- factor(Z, levels = 0:P)
    df <- data.frame(hash_matches)
    matches <- df %>%
      group_by(hash_matches, .drop = F) %>%
      count() %>%
      filter(hash_matches != 0) %>%
      pull()

    M.SAMPS[,s] <- m
    U.SAMPS[,s] <- u
    L.SAMPS[s] <- L
    PI.SAMPS[s] <- pi

    if(show_progress){
      if (s %% (S / 100) == 0) {
        flush.console()
        cat("\r", paste("Simulation", ": ", s / (S / 100), "% complete", sep = ""))
      }
    }
  }

  Z.SAMPS[Z.SAMPS == 0] <- n1 + 1

  list(Z = Z.SAMPS,
       m = M.SAMPS,
       u = U.SAMPS,
       overlap = L.SAMPS,
       pi = PI.SAMPS)

}

# gibbs_efficient <- function(hash, m_prior = 1, u_prior = 1,
#                        alpha = 1, beta = 1, S = 1000, burn = round(S * .1),
#                        show_progress = T){
#   # Implements bipartite record linkage with BK Sampling Mechanism
#   #
#   # Arguments
#   # comparisons = list calculated from from BRL::compareRecords
#   # m.prior = prior distribution for m parameters
#   # u.prior= prior distribution for u parameters
#   # alpha = first parameter of prior for linkage probability
#   # beta = second parameter of prior for linkage probability
#   # S = number of Gibbs iterations
#   # burn = number of iterations to be discarded as burn-in
#   # show_progress = set to false to show simulation progress
#
#   n1 <- hash$n1
#   n2 <- hash$n2
#   field_marker <- hash$field_marker
#
#   unique_patterns <- hash$ohe
#   pattern_counts <- hash$total_counts
#   P <- nrow(unique_patterns)
#   counts_by_rec <- hash$pattern_counts_by_record
#   hash_to_file_1 <-hash$hash_to_file_1
#
#   #candidates_P <- 1:(P+1)
#   candidates_P <- 0:P
#   Z.SAMPS <- matrix(NA, nrow = n2, ncol = S)
#   M.SAMPS <- matrix(NA, nrow = length(field_marker), ncol = S)
#   U.SAMPS <- matrix(NA, nrow = length(field_marker), ncol = S)
#   L.SAMPS <- vector(length = S)
#   PI.SAMPS <- vector(length = S)
#   Z <- rep(0, n2)
#   L <- 0
#
#   m <- u <- rep(0, length(field_marker))
#   matches <- rep(0,P)
#   #set.seed(1)
#
#   # Gibbs
#   for(s in 1:S){
#
#     AZ <- sweep(unique_patterns, MARGIN = 1, STAT = matches, FUN = "*") %>%
#       colSums() %>%
#       unname()
#
#     nonmatches <- pattern_counts - matches
#
#     BZ <- sweep(unique_patterns, MARGIN = 1, STAT = nonmatches, FUN = "*") %>%
#       colSums() %>%
#       unname()
#
#
#     m_post <- m_prior + AZ
#     u_post <- u_prior + BZ
#
#     m_post <- split(m_post, field_marker)
#     m <- as.vector(unlist(sapply(m_post, function(x){
#       prob <- MCMCpack::rdirichlet(1, x)
#       prob/sum(prob)
#     })))
#
#     u_post <- split(u_post, field_marker)
#     u <- as.vector(unlist(sapply(u_post, function(x){
#       prob <- MCMCpack::rdirichlet(1, x)
#       prob/sum(prob)
#     })))
#
#
#     ratio <- (log(m) - log(u)) %>%
#       rep(., P) %>%
#       matrix(., nrow = P, byrow = TRUE)
#
#     unique_weights <- exp(rowSums(ratio * unique_patterns, na.rm = TRUE))
#
#     hash_weights <- lapply(counts_by_rec, function(x){
#       x * unique_weights
#     })
#
#     pi <- rbeta(1, L + alpha, n2 - L + beta)
#     Z <- unname(sapply(hash_weights, function(x){
#       sample(candidates_P, 1, prob = c(1 - pi, x * pi / n1))
#     }))
#     #L <- sum(Z < P + 1)
#     L <- sum(Z > 0)
#     hash_matches <- factor(Z, levels = 0:P)
#     df <- data.frame(hash_matches)
#     matches <- df %>%
#       group_by(hash_matches, .drop = F) %>%
#       count() %>%
#       filter(hash_matches != 0) %>%
#       pull()
#
#     Z.SAMPS[,s] <- Z
#     M.SAMPS[,s] <- m
#     U.SAMPS[,s] <- u
#     L.SAMPS[s] <- L
#     PI.SAMPS[s] <- pi
#
#     if(show_progress){
#       if (s %% (S / 100) == 0) {
#         flush.console()
#         cat("\r", paste("Simulation", ": ", s / (S / 100), "% complete", sep = ""))
#       }
#     }
#   }
#
#   final_gibbs <- apply(Z.SAMPS, 2, function(z){
#     purrr::imap(z, ~ if(.x == 0) {
#       return(0)
#       } else {
#       sample_with_1(hash_to_file_1[[.y]][[.x]], 1)
#       }) %>%
#       unlist()
#   })
#
#   final_gibbs[final_gibbs == 0] <- n1 + 1
#
#   list(Z = final_gibbs,
#        m = M.SAMPS,
#        u = U.SAMPS,
#        overlap = L.SAMPS,
#        pi = PI.SAMPS)
#
# }


# gibbs_efficient_brian <- function(hash, m_prior = 1, u_prior = 1,
#                                   alpha = 1, beta = 1, S = 1000, burn = round(S * .1),
#                                   show_progress = T){
#   # Implements bipartite record linkage with BK Sampling Mechanism
#   #
#   # Arguments
#   # comparisons = list calculated from from BRL::compareRecords
#   # m.prior = prior distribution for m parameters
#   # u.prior= prior distribution for u parameters
#   # alpha = first parameter of prior for linkage probability
#   # beta = second parameter of prior for linkage probability
#   # S = number of Gibbs iterations
#   # burn = number of iterations to be discarded as burn-in
#   # show_progress = set to false to show simulation progress
#
#   n1 <- hash$n1
#   n2 <- hash$n2
#   field_marker <- hash$field_marker
#
#   unique_patterns <- hash$ohe
#   pattern_counts <- hash$total_counts
#   P <- nrow(unique_patterns)
#   counts_by_rec <- hash$pattern_counts_by_record
#   hash_to_file_1 <-hash$hash_to_file_1
#
#   #candidates_P <- 1:(P+1)
#   candidates_P <- 0:P
#   Z.SAMPS <- matrix(NA, nrow = n2, ncol = S)
#   M.SAMPS <- matrix(NA, nrow = length(field_marker), ncol = S)
#   U.SAMPS <- matrix(NA, nrow = length(field_marker), ncol = S)
#   L.SAMPS <- vector(length = S)
#   PI.SAMPS <- vector(length = S)
#   Z <- rep(0, n2)
#   L <- 0
#
#   m <- u <- rep(0, length(field_marker))
#   matches <- rep(0,P)
#   #set.seed(1)
#
#   # Gibbs
#   for(s in 1:S){
#
#     AZ <- sweep(unique_patterns, MARGIN = 1, STAT = matches, FUN = "*") %>%
#       colSums() %>%
#       unname()
#
#     nonmatches <- pattern_counts - matches
#
#     BZ <- sweep(unique_patterns, MARGIN = 1, STAT = nonmatches, FUN = "*") %>%
#       colSums() %>%
#       unname()
#
#
#     m_post <- m_prior + AZ
#     u_post <- u_prior + BZ
#
#     m_post <- split(m_post, field_marker)
#     m <- as.vector(unlist(sapply(m_post, function(x){
#       prob <- MCMCpack::rdirichlet(1, x)
#       prob/sum(prob)
#     })))
#
#     u_post <- split(u_post, field_marker)
#     u <- as.vector(unlist(sapply(u_post, function(x){
#       prob <- MCMCpack::rdirichlet(1, x)
#       prob/sum(prob)
#     })))
#
#
#     ratio <- (log(m) - log(u)) %>%
#       rep(., P) %>%
#       matrix(., nrow = P, byrow = TRUE)
#
#     unique_weights <- exp(rowSums(ratio * unique_patterns, na.rm = TRUE))
#
#     hash_weights <- lapply(counts_by_rec, function(x){
#       x * unique_weights
#     })
#
#     pi <- rbeta(1, L + alpha, n2 - L + beta)
#     Z <- unname(sapply(hash_weights, function(x){
#       sample(candidates_P, 1, prob = c(1 - pi, x * pi / n1))
#     }))
#     #L <- sum(Z < P + 1)
#     L <- sum(Z > 0)
#
#     records <- sapply(1:n2, function(j){
#       if(Z[j] == 0){
#         record = n1 + 1
#       } else {
#         index <- ceiling(runif(1) * counts_by_rec[[j]][Z[j]])
#         record <- hash_to_file_1[[j]][[Z[j]]][index]
#       }
#       record
#     })
#     # records <- purrr::imap(Z, ~ if(.x == 0) {
#     #   return(n1 + 1)
#     # } else {
#     #   #sample_with_1(hash_to_file_1[[.y]][[.x]], 1)
#     #   index <- ceiling(runif(1) * counts_by_rec[[.y]][.x])
#     #   record <- hash_to_file_1[[.y]][[.x]][index]
#     #   return(record)
#     # })
#
#     hash_matches <- factor(Z, levels = 0:P)
#     df <- data.frame(hash_matches)
#     matches <- df %>%
#       group_by(hash_matches, .drop = F) %>%
#       count() %>%
#       filter(hash_matches != 0) %>%
#       pull()
#
#     Z.SAMPS[,s] <- records
#     M.SAMPS[,s] <- m
#     U.SAMPS[,s] <- u
#     L.SAMPS[s] <- L
#     PI.SAMPS[s] <- pi
#
#     if(show_progress){
#       if (s %% (S / 100) == 0) {
#         flush.console()
#         cat("\r", paste("Simulation", ": ", s / (S / 100), "% complete", sep = ""))
#       }
#     }
#   }
#
#   list(Z = Z.SAMPS,
#        m = M.SAMPS,
#        u = U.SAMPS,
#        overlap = L.SAMPS,
#        pi = PI.SAMPS)
#
# }
