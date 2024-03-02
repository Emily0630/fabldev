
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

