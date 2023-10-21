#' @export
#'
vi_estimate_links_mm <- function(out, hash, resolve = TRUE,
                                lFNM=1, lFM1=1, lFM2=2, lR=Inf){

  # - positive losses
  C0 <- (lFNM > 0) & (lFM1 > 0) & (lFM2 > 0) & (lR > 0)
  # - conditions of Theorem 1 of Sadinle (2017)
  C1 <- (lR == Inf) & (lFNM <= lFM1) & (lFNM + lFM1 <= lFM2)
  # - conditions of Theorem 2 of Sadinle (2017)
  C2 <- ((lFM2 >= lFM1) & (lFM1 >= 2*lR)) | ((lFM1 >= lFNM) & (lFM2 >= lFM1 + lFNM))
  # - conditions of Theorem 3 of Sadinle (2017)
  C3 <- (lFM2 >= lFM1) & (lFM1 >= 2*lR) & (lFNM >= 2*lR)
  # check we can handle the specified losses
  if(!C0) stop("Losses need to be positive")
  if(!any(c(C1,C2,C3))) stop("Invalid configuration of losses")

  # tholdLink <- lFM1/(lFM1+lFNM) +
  #   (lFM2-lFM1-lFNM)*(1 - probNoLink - probMaxProbOption)/(lFM1+lFNM)
  tholdLink <- .5

  n2 <- hash$n2
  combination_counts <- hash$combination_counts
  pattern_weights <- out$pattern_weights
  C <- out$C
  pattern_finder <- hash$pattern_finder
  P <- length(pattern_finder)
  # combo_probs <- lapply(1:n2, function(j){
  #   out$pattern_weights/out$C[j]
  # })

  record <- 1

  matches <- lapply(1:n2, function(x){
    counts <- combination_counts[[x]] %>% unlist()
    probs <- pattern_weights * counts/ C[x]
    probs_combined <- sapply(1:P, function(y){
      sum(probs[pattern_finder[[y]]])
    })
    matching_patterns <- which(probs_combined > tholdLink)
    matching_probs <- probs_combined[matching_patterns]


    matching_records <- lapply(matching_patterns, function(y){
      hash$hash_to_file_1[[x]][[y]]
    }) %>%
      unlist()

    if(length(matching_records) == 0){
      return(data.frame(id_2 = x,
                        id_1 = NA,
                        prob = NA))
    } else {

    return(data.frame(id_2 = x,
               id_1 = matching_records,
               prob = matching_probs))
    }
  }) %>%
    do.call(rbind, .) %>%
    filter(!is.na(prob))

  Z_hat <- matches %>%
    select(id_1, id_2)

  prob <- matches$prob

  return(list(Z_hat = Z_hat,
              prob = prob))



  # maxProbOption <- max_prob$record
  # probMaxProbOption <- max_prob$prob
  # probNoLink <- out$b_pi/out$C
  # maxProbOptionIsLink <- maxProbOption > 0
  #
  # if(C1){# if not using reject option and conditions of Theorem 1
  #
  #   Z_hat <- rep(0, n2)
  #   tholdLink <- lFM1/(lFM1+lFNM) +
  #     (lFM2-lFM1-lFNM)*(1 - probNoLink - probMaxProbOption)/(lFM1+lFNM)
  #   Z_hat[maxProbOptionIsLink & (probMaxProbOption > tholdLink)] <-
  #     maxProbOption[maxProbOptionIsLink & (probMaxProbOption > tholdLink)]
  #
  # }else{# if using reject option
  #   if(C3){# if conditions of Theorem 3 are satisfied
  #
  #     Z_hat <- rep(-1,n2) # represents the reject option
  #     tholdLink <- 1 - lR/lFM1 + (lFM2-lFM1)*(1 - probNoLink - probMaxProbOption)/lFM1
  #     Z_hat[maxProbOptionIsLink & (probMaxProbOption > tholdLink) ] <-
  #       maxProbOption[maxProbOptionIsLink & (probMaxProbOption > tholdLink) ]
  #     noLinkDec <- probNoLink > 1-lR/lFNM
  #     #Z_hat[noLinkDec] <- ((n1+1):(n1+n2))[noLinkDec]
  #     Z_hat[noLinkDec] <- 0
  #
  #   }else{ # Theorem 2
  #
  #     # compute equation (6) in Sadinle (2017)
  #     # tableLabels[-n1-1,] <- t( lFM2*(t(1-tableLabels[-n1-1,])-tableLabels[n1+1,]) +
  #     #                             lFM1*tableLabels[n1+1,] )
  #     # tableLabels[n1+1,] <- lFNM*(1-tableLabels[n1+1,])
  #     # # find the options with the marginal minimal loss
  #     # lossMinLossOption <- apply(tableLabels, 2, min)
  #     # minLossOption <- apply(tableLabels, 2, which.min)
  #     # noLinkDec <- minLossOption == n1+1
  #     # minLossOption[noLinkDec] <- ((n1+1):(n1+n2))[noLinkDec]
  #     # Z_hat <- rep(-1,n2) # represents the reject option
  #     # Z_hat[lossMinLossOption < lR] <- minLossOption[lossMinLossOption < lR]
  #
  #   }
  # }
  #
  # if(resolve == T){
  # double_matches <- Z_hat[duplicated(Z_hat) & Z_hat > 0]
  # if (lR == Inf){
  #   to_resolve <- unlist(lapply(double_matches, function(x){
  #     dfB_options <- which(Z_hat == x)
  #     dfB_probs <- probMaxProbOption[dfB_options]
  #     non_matches <- dfB_options[-which.max(dfB_probs)]
  #     non_matches
  #   }))
  #   Z_hat[to_resolve] <- 0
  # } else {
  #   to_resolve <- unlist(lapply(double_matches, function(x){
  #     dfB_options <- which(Z_hat == x)
  #     dfB_options
  #   }))
  #   Z_hat[to_resolve] <- -1
  # }
  # }
  #
  # list(id_2 = Z_hat,
  #      prob = probMaxProbOption)

}

