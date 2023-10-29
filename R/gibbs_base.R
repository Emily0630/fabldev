#' @export
#'
gibbs_base <- function(comparisons, m_prior = 1, u_prior = 1,
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

  fields <- length(comparisons[[4]])
  n1 <- comparisons[[2]]; n2 <- comparisons[[3]]
  indicators_raw <-comparisons[[1]]

  field_marker <- as.vector(unlist(sapply(1:fields, function(x){
    rep(x, comparisons[[4]][x])
  })))

  ids <- expand.grid(1:n1, 1:n2)
  candidates <- 1:n1
  Z_samps <- matrix(NA, nrow = n2, ncol = S)
  m_samps <- matrix(NA, nrow = dim(indicators_raw)[2], ncol = S)
  u_samps <- matrix(NA, nrow = dim(indicators_raw)[2], ncol = S)
  L.SAMPS <- vector(length = S)
  Z.temp <- factor(rep(0, n1*n2), c(0,1))
  Z <- rep(n1+1, n2)
  L <- 0
  AZ <- BZ <- rep(0, length(field_marker))
  m <- u <- rep(0, length(field_marker))

  indicators <- data.frame(indicators_raw, Z.temp)

  # Gibbs
  for(s in 1:S){
    counts <- indicators %>%
      group_by(Z.temp, .drop = FALSE) %>%
      summarise(across(.cols = contains("X"),
                       .fns = sum))

    AZ<- counts %>%
      filter(Z.temp == 1) %>%
      select(contains("X")) %>%
      unlist(use.names = FALSE)

    BZ <- counts %>%
      filter(Z.temp == 0) %>%
      select(contains("X")) %>%
      unlist(use.names = FALSE)


    m.post <- m_prior + AZ
    u.post <- u_prior + BZ

    m.post <- split(m.post, field_marker)
    m <- as.vector(unlist(sapply(m.post, function(x){
      prob <- MCMCpack::rdirichlet(1, x)
      prob/sum(prob)
    })))

    u.post <- split(u.post, field_marker)
    u <- as.vector(unlist(sapply(u.post, function(x){
      prob <- MCMCpack::rdirichlet(1, x)
      prob/sum(prob)
    })))


    ratio <- (log(m) - log(u)) %>%
      rep(., n1 * n2) %>%
      matrix(., nrow = n1 *n2, byrow = TRUE)

    weights <- exp(rowSums(ratio * indicators_raw, na.rm = TRUE))

    weights <- split(weights, ids[,2])
    offset <- n1*(n2 - L + beta)/(L+ alpha)

    Z <- unname(sapply(weights, function(x){
      sample(c(candidates, n1 + 1), 1, prob = c(x, offset))
    }))

    Z.temp <- factor(as.vector(sapply(Z, function(x){
      if(x < n1 + 1){
        vec <- rep(0, n1)
        vec[x] <- 1
        vec
      }else{
        rep(0, n1)
      }
    })), c(0,1))

    L <- sum(Z < n1 + 1)
    indicators$Z.temp <- Z.temp

    Z_samps[,s] <- Z
    m_samps[,s] <- m
    u_samps[,s] <- u
    L.SAMPS[s] <- L

    if(show_progress){
      if (s %% (S / 100) == 0) {
        flush.console()
        cat("\r", paste("Simulation", ": ", s / (S / 100), "% complete", sep = ""))
      }
    }
  }
  Z_samps <- Z_samps[,-(1:burn)]
  m_samps <- m_samps[,-(1:burn)]
  u_samps <- u_samps[,-(1:burn)]

  list(Z = Z_samps,
       m = m_samps,
       u = u_samps)

}


