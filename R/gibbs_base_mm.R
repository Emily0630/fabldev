#' @export
#'
gibbs_base_mm <- function(comparisons, m_prior = 1, u_prior = 1,
                          alpha = 1, beta = 1, S = 1000, burn = 100,
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

  nonmatch_df <- data.frame(id_1 = 0,
                            id_2 = 1:n2,
                            rn = NA)

  ids_df <- expand.grid(1:n1, 1:n2) %>%
    setNames(., c("id_1", "id_2")) %>%
    mutate(rn = row_number())

  candidates <- 1:n1
  Z_list <- vector("list", S)
  m_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
  u_samps <- matrix(NA, nrow = length(field_marker), ncol = S)

  pi_samps <- list()
  n_sampled_list <- list()
  Z.temp <- factor(rep(0, n1*n2), c(0,1))

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
    #weights <- split(weights, ids[, 2])

    #ids_df$weight[ids_df$id_1 > 0] <- weights
    ids_df$weight <- weights
    #rn_long <- ids_df$rn

    k  = 1
    n_sampled_vec <- c()
    pi_vec <- c()
    n_possible_vec <- n2
    #n_possible <- n2
    matchable <- 1:n2
    removed_set <- c()
    available <- ids_df


    while(TRUE){
      if(s == 1){
        n_prior = 0
      } else if(length(n_sampled_list[[s - 1]]) < k){
        n_prior = 0
      } else {
        n_prior <- n_sampled_list[[s - 1]][k]
      }

      n_possible <- n_possible_vec[k]



      pi <- rbeta(1, n_prior + alpha, n_possible - n_prior + beta)
      offset <- (n1 - k + 1) * (1 - pi) / pi
      #ids_df$weight[ids_df$id_1 == 0] <- offset
      available_split <- available %>%
        group_split(id_2)


      # sampled <- ids_df %>%
      #   filter(id_2 %in% matchable,
      #          !(rn %in% removed_set)) %>%
      #   group_by(id_2) %>%
      #   sample_n(1, weight = weight) %>%
      #   ungroup()

      sampled_rows <- sapply(matchable, function(x){
        sample(c(0, available_split[[x]]$rn),
               1,
               prob = c(offset, available_split[[x]]$weight))
      }) %>%
        unname()

      removed <- sampled_rows[sampled_rows > 0]
      matchable <- ids_df %>%
        filter(rn %in% removed) %>%
        select(id_2) %>%
        pull()



      removed_set <- c(removed_set, removed)
      n_possible_vec <- c(n_possible_vec, length(matchable))

      available <- available %>%
        filter(!(rn %in% removed_set))

      n_sampled_vec <- c(n_sampled_vec, length(matchable))
      pi_vec <- c(pi_vec, pi)

      k = k + 1

      if(length(matchable) == 0){
        break
      }

    }

    n_sampled_list[[s]] <- n_possible_vec[-1]
    pi_samps[[s]] <- pi_vec

    # Z_compact[[s]] <- ids_df %>%
    #   filter(rn %in% removed_set) %>%
    #   tidyr::complete(id_2 = 1:n2) %>%
    #   select(id_1, id_2) %>%
    #   nest_by(id_2, .key = "id_1") %>%
    #   ungroup() %>%
    #   select(id_1)

    Z <- ids_df %>%
      filter(rn %in% removed_set) %>%
      select(id_1, id_2)

    Z.temp <- rep(0, (n1 * n2))
    Z.temp[removed_set] <- 1



    indicators$Z.temp <- Z.temp

    Z_list[[s]] <- Z
    m_samps[,s] <- m
    u_samps[,s] <- u

    if(show_progress){
      if (s %% (S / 100) == 0) {
        flush.console()
        cat("\r", paste("Simulation", ": ", s / (S / 100), "% complete", sep = ""))
      }
    }
  }

  # Z_samps <- lapply(Z_list, function(y){
  #   y %>%
  #     tidyr::complete(id_2 = 1:n2) %>%
  #     nest_by(id_2, .key = "id_1") %>%
  #     ungroup() %>%
  #     select(id_1)
  # }) %>%
  #   do.call(cbind, .)


  #Z_samps <- Z_samps[,-(1:burn)]

  Z_samps <- lapply((burn + 1):S, function(x){
    Z_list[[x]]
  })

  m_samps <- m_samps[,-(1:burn)]
  u_samps <- u_samps[,-(1:burn)]

  list(Z = Z_samps,
       m = m_samps,
       u = u_samps)

}




# gibbs_base_mm <- function(comparisons, m_prior = 1, u_prior = 1,
#                           alpha = 1, beta = 1, S = 1000, burn = 100,
#                           show_progress = T){
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
#   fields <- length(comparisons[[4]])
#   levels <-
#     n1 <- comparisons[[2]]; n2 <- comparisons[[3]]
#   indicators_raw <-comparisons[[1]]
#
#   field_marker <- as.vector(unlist(sapply(1:fields, function(x){
#     rep(x, comparisons[[4]][x])
#   })))
#
#   ids <- expand.grid(1:n1, 1:n2)
#
#   nonmatch_df <- data.frame(id_1 = 0,
#                             id_2 = 1:n2,
#                             rn = NA)
#
#   ids_df <- expand.grid(1:n1, 1:n2) %>%
#     setNames(., c("id_1", "id_2")) %>%
#     mutate(rn = row_number()) %>%
#     bind_rows(nonmatch_df)
#
#
#
#   candidates <- 1:n1
#   #Z_samps <- matrix(NA, nrow = n2, ncol = S)
#   Z_samps <- matrix(NA, nrow = n1*n2, ncol = S)
#   Z_compact <- list()
#   m_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
#   u_samps <- matrix(NA, nrow = length(field_marker), ncol = S)
#   L.SAMPS <- vector(length = S)
#   pi_samps <- list()
#   n_sampled_list <- list()
#   Z.temp <- factor(rep(0, n1*n2), c(0,1))
#   Z <- rep(n1+1, n2)
#   L <- 0
#   AZ <- BZ <- rep(0, length(field_marker))
#   m <- u <- rep(0, length(field_marker))
#
#   indicators <- data.frame(indicators_raw, Z.temp)
#
#   # Gibbs
#   for(s in 1:S){
#     counts <- indicators %>%
#       group_by(Z.temp, .drop = FALSE) %>%
#       summarise(across(.cols = contains("X"),
#                        .fns = sum))
#
#     AZ<- counts %>%
#       filter(Z.temp == 1) %>%
#       select(contains("X")) %>%
#       unlist(use.names = FALSE)
#
#     BZ <- counts %>%
#       filter(Z.temp == 0) %>%
#       select(contains("X")) %>%
#       unlist(use.names = FALSE)
#
#
#     m.post <- m_prior + AZ
#     u.post <- u_prior + BZ
#
#     m.post <- split(m.post, field_marker)
#     m <- as.vector(unlist(sapply(m.post, function(x){
#       prob <- MCMCpack::rdirichlet(1, x)
#       prob/sum(prob)
#     })))
#
#     u.post <- split(u.post, field_marker)
#     u <- as.vector(unlist(sapply(u.post, function(x){
#       prob <- MCMCpack::rdirichlet(1, x)
#       prob/sum(prob)
#     })))
#
#     ratio <- (log(m) - log(u)) %>%
#       rep(., n1 * n2) %>%
#       matrix(., nrow = n1 *n2, byrow = TRUE)
#
#     weights <- exp(rowSums(ratio * indicators_raw, na.rm = TRUE))
#     #weights <- split(weights, ids[, 2])
#
#     ids_df$weight[ids_df$id_1 > 0] <- weights
#
#     k  = 1
#     n_sampled_vec <- c()
#     pi_vec <- c()
#     n_possible_vec <- n2
#     #n_possible <- n2
#     matchable <- 1:n2
#     removed_set <- c()
#
#
#     while(TRUE){
#       if(s == 1){
#         n_prior = 0
#       } else if(length(n_sampled_list[[s - 1]]) < k){
#         n_prior = 0
#       } else {
#         n_prior <- n_sampled_list[[s - 1]][k]
#       }
#
#       n_possible <- n_possible_vec[k]
#
#
#
#       pi <- rbeta(1, n_prior + alpha, n_possible - n_prior + beta)
#       offset <- (n1 - k + 1) * (1 - pi) / pi
#       ids_df$weight[ids_df$id_1 == 0] <- offset
#       #
#       # weights_df <- ids %>%
#       #   filter(!(row_number %in% removed)) %>%
#       #   nest_by(id_2)
#
#       sampled <- ids_df %>%
#         filter(id_2 %in% matchable,
#                !(rn %in% removed_set)) %>%
#         group_by(id_2) %>%
#         sample_n(1, weight = weight) %>%
#         ungroup()
#
#       #   sampled_rows <- sapply(matchable, function(x){
#       #     sample(c(0, rn_split[[x]][available[[x]]]),
#       #            1,
#       #            prob = c(offset, weights[[x]][available[[x]]]))
#       # }) %>%
#       #     unname()
#
#       matchable <- sampled %>%
#         filter(!is.na(rn)) %>%
#         select(id_2) %>%
#         pull()
#       removed <- sampled %>%
#         filter(!is.na(rn)) %>%
#         select(rn) %>%
#         pull()
#       removed_set <- c(removed_set, removed)
#       n_possible_vec <- c(n_possible_vec, length(matchable))
#
#       n_sampled_vec <- c(n_sampled_vec, length(matchable))
#       pi_vec <- c(pi_vec, pi)
#
#       k = k + 1
#
#       if(length(matchable) == 0){
#         break
#       }
#
#     }
#
#     n_sampled_list[[s]] <- n_possible_vec[-1]
#     pi_samps[[s]] <- pi_vec
#
#     Z_compact[[s]] <- ids_df %>%
#       filter(rn %in% removed_set) %>%
#       tidyr::complete(id_2 = 1:n2) %>%
#       select(id_1, id_2) %>%
#       nest_by(id_2, .key = "id_1") %>%
#       ungroup() %>%
#       select(id_1)
#
#     Z.temp <- rep(0, (n1 * n2))
#     Z.temp[removed_set] <- 1
#
#
#     # Z.temp <- factor(as.vector(sapply(Z, function(x){
#     #   if(x < n1 + 1){
#     #     vec <- rep(0, n1)
#     #     vec[x] <- 1
#     #     vec
#     #   }else{
#     #     rep(0, n1)
#     #   }
#     # })), c(0,1))
#
#     # L <- sum(Z > 0)
#     indicators$Z.temp <- Z.temp
#
#     #Z_samps[,s] <- Z
#     Z_samps[,s] <- Z.temp
#     m_samps[,s] <- m
#     u_samps[,s] <- u
#     #L.SAMPS[s] <- L
#
#     if(show_progress){
#       if (s %% (S / 100) == 0) {
#         flush.console()
#         cat("\r", paste("Simulation", ": ", s / (S / 100), "% complete", sep = ""))
#       }
#     }
#   }
#
#   Z_compact <- do.call(cbind, Z_compact)
#   Z_compact <- Z_compact[,-(1:burn)]
#   Z_samps <- Z_samps[,-(1:burn)]
#   #L.SAMPS <- L.SAMPS[-(1:burn)]
#   m_samps <- m_samps[,-(1:burn)]
#   u_samps <- u_samps[,-(1:burn)]
#
#   list(Z = Z_samps,
#        Z_compact = Z_compact,
#        m = m_samps,
#        u = u_samps)
#
# }



