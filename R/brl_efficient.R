#' @export
#'
brl_efficient <- function(hash, rejection = F, m_prior = 1, u_prior = 1,
                          alpha = 1, beta = 1, S = 1000, burn = 100,
                          show_progress = T, seed = 0){
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

  set.seed(seed)

  n1 <- hash$n1
  n2 <- hash$n2
  field_marker <- hash$field_marker

  unique_patterns <- hash$ohe
  pattern_counts <- hash$total_counts
  P <- nrow(unique_patterns)
  counts_by_rec <- hash$pattern_counts_by_record
  hash_to_file_1 <-hash$hash_to_file_1

  #candidates_P <- 1:(P+1)
  candidates_P <- 0:P
  Z.SAMPS <- matrix(NA, nrow = n2, ncol = S)
  M.SAMPS <- matrix(NA, nrow = length(field_marker), ncol = S)
  U.SAMPS <- matrix(NA, nrow = length(field_marker), ncol = S)
  L.SAMPS <- vector(length = S)
  Z <- rep(0, n2)
  Z_pattern <- rep(0, n2)
  L <- 0

  m <- u <- rep(0, length(field_marker))
  matches <- rep(0, P)
  #already_matched <- rep(0, n1)
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


    for(j in sample(1:n2)){

      if(Z[j] > 0){
        L <- L - 1
      }
      Z[j] <- 0

      if(!rejection){
        already_matched <- Z[Z > 0]
        temp_hash <- lapply(hash_to_file_1[[j]], function(x){
          x[!(x %in% already_matched)]
        })
        temp_counts <- sapply(temp_hash, length)
        temp_weights <- unique_weights * temp_counts

        Z_pattern[j] <-
          sample(candidates_P, 1,
                 prob = c((n1 - L) * (n2 - L - 1 + beta) / (L + alpha),
                          temp_weights))

        Z_record <- if(Z_pattern[j] == 0){
          0
        } else {
          sample_with_1(temp_hash[[Z_pattern[j]]], 1)
        }
        Z[j] <- Z_record
      }
      else{
        hash_weights <- counts_by_rec[[j]] * unique_weights
        probs <- c((n1 - L) * (n2 - L - 1 + beta) / (L + alpha),
                   hash_weights)
        flag <- 1
        while(flag == 1){
          pattern <- sample(candidates_P, 1, prob = probs)
          if(pattern == 0){
            Z[j] <- 0
            flag <- 0
          }
          else{
            index <- ceiling(runif(1) * counts_by_rec[[j]][pattern])
            i <- hash_to_file_1[[j]][[pattern]][index]
            if(!(i %in% Z)){
              Z[j] <- i
              flag <- 0
            }
          }
        }
      }

      if(Z[j] > 0){
        L <- L + 1
      }
    }
    hash_matches <- factor(Z_pattern, levels = 0:P)
    df <- data.frame(hash_matches)
    matches <- df %>%
      group_by(hash_matches, .drop = F) %>%
      count() %>%
      filter(hash_matches != 0) %>%
      pull()

    Z.SAMPS[,s] <- Z
    M.SAMPS[,s] <- m
    U.SAMPS[,s] <- u
    L.SAMPS[s] <- L

    if(show_progress){
      if (s %% (S / 100) == 0) {
        flush.console()
        cat("\r", paste("Simulation", ": ", s / (S / 100), "% complete",
                        sep = ""))
      }
    }
  }

  # final_gibbs <- apply(Z.SAMPS, 2, function(z){
  #   purrr::imap(z, ~ if(.x == 0) {
  #     return(0)
  #     } else {
  #     sample_with_1(hash_to_file_1[[.y]][[.x]], 1)
  #     }) %>%
  #     unlist()
  # })

  #final_gibbs[final_gibbs == 0] <- n1 + 1

  Z.SAMPS[Z.SAMPS == 0] <- n1 + 1

  list(Z = Z.SAMPS,
       m = M.SAMPS,
       u = U.SAMPS,
       overlap = L.SAMPS)

}
