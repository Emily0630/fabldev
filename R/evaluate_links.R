#' @export
#'
evaluate_links <- function(Z_hat, Z_true, n1){
  # Z_hat = Bayes Estimate of linkage structure (from BRL)
  # Z_true = True linkage structure
  # n1 = size of larger file

  nLinks <- sum(Z_hat <= n1 & Z_hat > 0)
  nMatches <- sum(Z_true <= n1 & Z_true > 0)
  nCorrectLinks <- sum(Z_hat[Z_hat<=n1 & Z_hat > 0]==Z_true[Z_hat<=n1 & Z_hat > 0])
  recall <- nCorrectLinks/nMatches
  precision <- nCorrectLinks/nLinks
  fmeasure <- 2 * (recall * precision) / (recall + precision)
  eval <- c(recall, precision, fmeasure)
  names(eval) <- c("Recall", "Precision", "Fmeasure")
  eval
}

# evaluate_links_mm <- function(Z_hat, Z_true, n1){
#   # Z_hat = Bayes Estimate of linkage structure (from BRL)
#   # Z_true = True linkage structure
#   # n1 = size of larger file
#
#   n2 <- nrow(Z_true)
#   n_nonlinks <- Z_hat %>%
#     tidyr::unnest(cols = c(data)) %>%
#     filter(is.na(id_1)) %>%
#     nrow()
#   n_links <- n2 - n_nonlinks
#
#   n_nonmatches <- Z_true %>%
#     tidyr::unnest(cols = c(data)) %>%
#     filter(is.na(id_1)) %>%
#     nrow()
#
#   n_matches <- n2 - n_nonmatches
#   n_correct_links <- sapply(1:n2, function(x){
#     vec1 = unlist(Z_true[x, ]$data)
#     vec2 = unlist(Z_hat[x, ]$data)
#     if(length(vec1) != length(vec2)){
#       FALSE
#     } else{
#       all(vec1 == vec2)
#     }
#   }) %>%
#     sum(na.rm = T)
#
#   recall <- n_correct_links / n_matches
#   precision <- n_correct_links/n_links
#   fmeasure <- 2 * (recall * precision) / (recall + precision)
#   eval <- c(recall, precision, fmeasure)
#   names(eval) <- c("Recall", "Precision", "Fmeasure")
#   eval
# }

# TODO: Evalution of single match vs multiple match truth
# TODO: Pair based evaluation metric

evaluate_links_mm <- function(Z_hat, Z_true, n1){
  # Z_hat = Bayes Estimate of linkage structure (from BRL)
  # Z_true = True linkage structure
  # n1 = size of larger file

  if(typeof(Z_hat) == "double"){
    Z_hat <- data.frame(id_1 = Z_hat,
                        id_2 = 1:n2) %>%
      filter(id_1 > 0)

  } else {

  }

  n2 <- unique(Z_true$id_2) %>%
    length()
  # n_nonlinks <- Z_hat %>%
  #   tidyr::unnest(cols = c(data)) %>%
  #   filter(is.na(id_1)) %>%
  #   nrow()
  n_links <- nrow(Z_hat)

  n_matches <- Z_true %>%
    filter(!is.na(id_1)) %>%
    nrow(.)

  Z_hat_char <- Z_hat %>%
    tidyr::unite(col = "pair", sep = ",") %>%
    pull()

  Z_true_char <- Z_true %>%
    tidyr::unite(col = "pair", sep = ",") %>%
    pull()

  n_correct_links <- sapply(Z_hat_char, function(x){
    x %in% Z_true_char
  }) %>%
    sum()

  recall <- n_correct_links / n_matches
  precision <- n_correct_links/n_links
  fmeasure <- 2 * (recall * precision) / (recall + precision)
  eval <- c(recall, precision, fmeasure)
  names(eval) <- c("Recall", "Precision", "Fmeasure")
  eval
}

