estimate_links_mm <- function(Z_samps, n1, l_FNM=1, l_FM1=1, l_FM2=2, l_R=Inf,
                              resolve = T){
  #
  # This is a complete copy of "linkrecords" from BRL, only modified
  # so that it passes on the posterior link probabilities
  #
  #
  #

  # control the input
  #if(!is.matrix(Z_samps)) stop("Z_samps should be a matrix")
  #n2 <- nrow(Z_samps)
  # make sure the labels in Z_samps are within the expected range
  #if(max(Z_samps) > n1 + n2) stop("Labels in Z_samps exceed n1+n2")
  # - positive losses
  C0 <- (l_FNM > 0) & (l_FM1 > 0) & (l_FM2 > 0) & (l_R > 0)
  # - conditions of Theorem 1 of Sadinle (2017)
  C1 <- (l_R == Inf) & (l_FNM <= l_FM1) & (l_FNM + l_FM1 <= l_FM2)
  # - conditions of Theorem 2 of Sadinle (2017)
  C2 <- ((l_FM2 >= l_FM1) & (l_FM1 >= 2*l_R)) | ((l_FM1 >= l_FNM) & (l_FM2 >= l_FM1 + l_FNM))
  # - conditions of Theorem 3 of Sadinle (2017)
  C3 <- (l_FM2 >= l_FM1) & (l_FM1 >= 2*l_R) & (l_FNM >= 2*l_R)
  # check we can handle the specified losses
  if(!C0) stop("Losses need to be positive")
  if(!any(c(C1,C2,C3))) stop("Invalid configuration of losses")

  #threshold <- l_FM1/(l_FM1+l_FNM)
  threshold <- 1/2
  temp <- Z_samps %>%
    do.call(rbind, .) %>%
    group_by(id_1, id_2) %>%
    count() %>%
    mutate(prob = n / length(Z_samps)) %>%
    filter(prob > threshold) %>%
    ungroup()

  Z_hat <- temp %>%
    select(id_1, id_2)

  prob <- temp %>%
    select(prob)

  double_matches <- Z_hat$id_1[duplicated(Z_hat$id_1)]

  if(resolve == TRUE & length(double_matches) > 0){
    if (l_R == Inf){
      to_resolve <- unlist(lapply(double_matches, function(x){
        df2_options <- which(Z_hat$id_1 == x)
        df2_probs <- probs[df2_options]
        non_matches <- df2_options[-which.max(df2_probs)]
        non_matches
      }))
      Z_hat <- Z_hat[-to_resolve, ]
      prob <- prob[-to_resolve, ]
    }
  }

  return(list(Z_hat = Z_hat,
              prob = prob))
}




