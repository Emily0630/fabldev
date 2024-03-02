#' @export
#'
estimate_links<- function(Z_samps, n1, l_FNM=1, l_FM1=1, l_FM2=2, l_R=Inf,
                          nonmatch_label = "zero", resolve = T){
  #
  # Adapted from BRL::linkrecords. See https://CRAN.R-project.org/package=BRL
  #

  n2 <- nrow(Z_samps)

  # - positive losses
  positive_losses <- (l_FNM > 0) & (l_FM1 > 0) & (l_FM2 > 0) & (l_R > 0)
  # - conditions of Theorem 1 of Sadinle (2017)
  C1 <- (l_R == Inf) & (l_FNM <= l_FM1) & (l_FNM + l_FM1 <= l_FM2)
  # - conditions of Theorem 2 of Sadinle (2017)
  theorem_2 <- ((l_FM2 >= l_FM1) & (l_FM1 >= 2*l_R)) | ((l_FM1 >= l_FNM) & (l_FM2 >= l_FM1 + l_FNM))
  # - conditions of Theorem 3 of Sadinle (2017)
  C3 <- (l_FM2 >= l_FM1) & (l_FM1 >= 2*l_R) & (l_FNM >= 2*l_R)
  # check we can handle the specified losses
  if(!positive_losses) stop("Losses need to be positive")
  if(!theorem_2) stop("Not yet configured for Sadinle (2017) Theorem 2")


  # temporarily replace all nonlink labels by n1+1
  Z_samps[Z_samps > n1+1] <- n1+1

    samps <- ncol(Z_samps)
    probs <- apply(Z_samps, 1, function(x){
      table(x)/samps
    })
    prob_no_link <- sapply(probs, function(x){
      1 - sum(x[names(x) != n1 + 1])
    })
    Z_hat <- rep(0, n2)
    best_match <- sapply(probs, function(x){
      names(which.max(x))
    }) %>%
      as.numeric()
    prob_best_match <- sapply(probs, function(x){
      max(x)
    })
    link_indicator <- best_match < n1 + 1

  if(l_R == Inf){# if not using reject option and conditions of Theorem 1

    if (nonmatch_label == "n_1 + j"){
      Z_hat <- (n1+1):(n1+n2)
    }
    if (nonmatch_label == "zero"){
      Z_hat <- rep(0, n2)
    }
    threshold <- l_FM1/(l_FM1+l_FNM) +
      (l_FM2-l_FM1-l_FNM)*(1 - prob_no_link - prob_best_match)/(l_FM1+l_FNM)

    Z_hat[link_indicator & (prob_best_match > threshold)] <-
      best_match[link_indicator & (prob_best_match > threshold)]

  }

  if(l_R < Inf){
    if(!theorem_2){
      Z_hat <- rep(-1, n2) # represents the reject option
      threshold <- 1 - l_R/l_FM1 +
        (l_FM2-l_FM1)*(1 - prob_no_link - prob_best_match) / l_FM1

      Z_hat[link_indicator & (prob_best_match > threshold) ] <-
        best_match[link_indicator & (prob_best_match > threshold) ]
      noLinkDec <- prob_no_link > 1 - l_R / l_FNM

      if (nonmatch_label == "n_1 + j"){
        Z_hat[noLinkDec] <- ((n1+1):(n1+n2))[noLinkDec]
      }
      if (nonmatch_label == "zero"){
        Z_hat[noLinkDec] <- 0
      }
    } else {
    # TODO: Write code for Theorem 2
    }
  }


  # Enforce one-to-one matching
  if (resolve == T){
    double_matches <- Z_hat[duplicated(Z_hat) & Z_hat > 0]
    if (l_R == Inf){
      to_resolve <- unlist(lapply(double_matches, function(x){
        dfB_options <- which(Z_hat == x)
        dfB_probs <- prob_best_match[dfB_options]
        non_matches <- dfB_options[-which.max(dfB_probs)]
        non_matches
      }))
      Z_hat[to_resolve] <- 0
    } else {
      to_resolve <- unlist(lapply(double_matches, function(x){
        dfB_options <- which(Z_hat == x)
        dfB_options
      }))
      Z_hat[to_resolve] <- -1
    }
  }

  return(list(Z_hat = Z_hat,
              prob = prob_best_match))
}
