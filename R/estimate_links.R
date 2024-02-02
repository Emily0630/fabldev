#' @export
#'
estimate_links<- function(Z_samps, n1, lFNM=1, lFM1=1, lFM2=2, lR=Inf,
                          nonmatch = "BK", matching = "single", resolve = T,
                          edition = "new"){
  #
  # Adapted from BRL::linkrecords. See https://CRAN.R-project.org/package=BRL
  #
  #
  #

  # control the input
  if(!is.matrix(Z_samps)) stop("Z_samps should be a matrix")
  n2 <- nrow(Z_samps)
  # make sure the labels in Z_samps are within the expected range
  if(max(Z_samps) > n1 + n2) stop("Labels in Z_samps exceed n1+n2")
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

  # temporarily replace all nonlink labels by n1+1
  Z_samps[Z_samps > n1+1] <- n1+1
  if(edition == "old"){
  tableLabels <- apply(Z_samps, 1, tabulate, nbins=max(Z_samps))
  tableLabels <- tableLabels/ncol(Z_samps)
  probNoLink <- tableLabels[n1+1,]
  # find marginal best option for each record based only on probability
  maxProbOption <- apply(tableLabels, 2, which.max)
  maxProbOption[maxProbOption==n1+1] <- (n1+1:n2)[maxProbOption==n1+1]
  probMaxProbOption <- apply(tableLabels, 2, max)
  maxProbOptionIsLink <- maxProbOption <= n1
  }

  if(edition == "new"){
    samps <- ncol(Z_samps)
    probs <- apply(Z_samps, 1, function(x){
      table(x)/samps
    })
    probNoLink <- sapply(probs, function(x){
      1 - sum(x[names(x) != n1 + 1])
    })
    Z_hat <- rep(0, n2)
    maxProbOption <- sapply(probs, function(x){
      names(which.max(x))
    }) %>%
      as.numeric()
    probMaxProbOption <- sapply(probs, function(x){
      max(x)
    })
    maxProbOptionIsLink <- maxProbOption < n1 + 1
    }

  if(C1){# if not using reject option and conditions of Theorem 1

    if (nonmatch == "Sadinle"){
      Z_hat <- (n1+1):(n1+n2)
    }
    if (nonmatch == "BK"){
      Z_hat <- rep(0, n2)
    }
    tholdLink <- lFM1/(lFM1+lFNM) +
      (lFM2-lFM1-lFNM)*(1 - probNoLink - probMaxProbOption)/(lFM1+lFNM)
    Z_hat[maxProbOptionIsLink & (probMaxProbOption > tholdLink)] <-
      maxProbOption[maxProbOptionIsLink & (probMaxProbOption > tholdLink)]

  }else{# if using reject option
    if(C3){# if conditions of Theorem 3 are satisfied

      Z_hat <- rep(-1,n2) # represents the reject option
      tholdLink <- 1 - lR/lFM1 + (lFM2-lFM1)*(1 - probNoLink - probMaxProbOption)/lFM1
      Z_hat[maxProbOptionIsLink & (probMaxProbOption > tholdLink) ] <-
        maxProbOption[maxProbOptionIsLink & (probMaxProbOption > tholdLink) ]
      noLinkDec <- probNoLink > 1-lR/lFNM

      if (nonmatch == "Sadinle"){
        Z_hat[noLinkDec] <- ((n1+1):(n1+n2))[noLinkDec]
      }
      if (nonmatch == "BK"){
        Z_hat[noLinkDec] <- 0
      }


    }else{ # Theorem 2

      # compute equation (6) in Sadinle (2017)
      tableLabels[-n1-1,] <- t( lFM2*(t(1-tableLabels[-n1-1,])-tableLabels[n1+1,]) +
                                  lFM1*tableLabels[n1+1,] )
      tableLabels[n1+1,] <- lFNM*(1-tableLabels[n1+1,])
      # find the options with the marginal minimal loss
      lossMinLossOption <- apply(tableLabels, 2, min)
      minLossOption <- apply(tableLabels, 2, which.min)
      noLinkDec <- minLossOption == n1+1
      minLossOption[noLinkDec] <- ((n1+1):(n1+n2))[noLinkDec]
      Z_hat <- rep(-1,n2) # represents the reject option
      Z_hat[lossMinLossOption < lR] <- minLossOption[lossMinLossOption < lR]

    }
  }

  if (resolve == T){

    double_matches <- Z_hat[duplicated(Z_hat) & Z_hat > 0]
    if (lR == Inf){
      to_resolve <- unlist(lapply(double_matches, function(x){
        dfB_options <- which(Z_hat == x)
        dfB_probs <- probMaxProbOption[dfB_options]
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
              prob = probMaxProbOption))
}

# estimate_links_mm<- function(Z_samps, n1, lFNM=1, lFM1=1, lFM2=2, lR=Inf,
#                           nonmatch = "BK"){
#   #
#   # This is a complete copy of "linkrecords" from BRL, only modified
#   # so that it passes on the posterior link probabilities
#   #
#   #
#   #
#
#   # control the input
#   #if(!is.matrix(Z_samps)) stop("Z_samps should be a matrix")
#   n2 <- nrow(Z_samps)
#   # make sure the labels in Z_samps are within the expected range
#   #if(max(Z_samps) > n1 + n2) stop("Labels in Z_samps exceed n1+n2")
#   # - positive losses
#   C0 <- (lFNM > 0) & (lFM1 > 0) & (lFM2 > 0) & (lR > 0)
#   # - conditions of Theorem 1 of Sadinle (2017)
#   C1 <- (lR == Inf) & (lFNM <= lFM1) & (lFNM + lFM1 <= lFM2)
#   # - conditions of Theorem 2 of Sadinle (2017)
#   C2 <- ((lFM2 >= lFM1) & (lFM1 >= 2*lR)) | ((lFM1 >= lFNM) & (lFM2 >= lFM1 + lFNM))
#   # - conditions of Theorem 3 of Sadinle (2017)
#   C3 <- (lFM2 >= lFM1) & (lFM1 >= 2*lR) & (lFNM >= 2*lR)
#   # check we can handle the specified losses
#   if(!C0) stop("Losses need to be positive")
#   if(!any(c(C1,C2,C3))) stop("Invalid configuration of losses")
#
#   #tholdLink <- lFM1/(lFM1+lFNM)
#   threshold <- 1/2
#
#   Z_hat <- lapply(1:n2, function(x){
#     Z_samps[x, ] %>%
#       unlist() %>%
#       unname() %>%
#       data.frame(id_1 = .) %>%
#       group_by(id_1) %>%
#       count() %>%
#       mutate(prob = n/ ncol(Z_samps),
#              id_2 = x) %>%
#       filter(prob > threshold) %>%
#       ungroup() %>%
#       select(-n, -prob)
#
#     # %>%
#     #   nest_by(id_2)
#   }) %>%
#     do.call(rbind, .) %>%
#     nest_by(id_2, .keep = F)
#
#   # Z_hat <- lapply(1:n2, function(x){
#   #   Z_samps[x, ] %>%
#   #     unlist() %>%
#   #     unname() %>%
#   #     data.frame(id_1 = .) %>%
#   #     group_by(id_1) %>%
#   #     count() %>%
#   #     mutate(id_2 = x,
#   #            prob = n/ ncol(Z_samps)) %>%
#   #     filter(prob > threshold) %>%
#   #     select(-n) %>%
#   #     ungroup()
#   # }) %>%
#   #   do.call(rbind, .)
#
#
#   return(Z_hat)
# }


estimate_links_mm <- function(Z_samps, n1, lFNM=1, lFM1=1, lFM2=2, lR=Inf,
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

  #threshold <- lFM1/(lFM1+lFNM)
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
    if (lR == Inf){
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




