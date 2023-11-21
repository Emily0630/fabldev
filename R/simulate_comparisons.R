#' Simulate comparison Vectors
#'
#' @export

simulate_comparisons <- function(m, u, levels, n1, n2, overlap,
                                 previous_matches = 0){
  field_marker <- unlist(lapply(1:length(levels), function(x){
    rep(x, levels[x])
  }))

  N <- n1 * n2
  ids <- expand.grid(1:n1, 1:n2)
  indicators <- matrix(NA, nrow = N, ncol = length(levels))

  df2matches <- seq_len(overlap)
  df1matches <- df2matches + previous_matches

  Ztrue <- rep(n1 + 1, n2)
  Ztrue[df2matches] <- df1matches

  match_index <- which(ids[,1]  == (ids[,2]+ previous_matches))[seq_len(overlap)]

  m_list <- split(m, field_marker)
  u_list <- split(u, field_marker)

  gamma_match <- sapply(m_list, function(x){
    sample(seq_along(x) - 1, overlap, replace = T, x)
  })

  gamma_nonmatch <- sapply(u_list, function(x){
    sample(seq_along(x) - 1, N - overlap, replace = T, x)
  })

  if(overlap == 0){
    indicators <- gamma_nonmatch
  } else {
  indicators[match_index,] <- gamma_match
  indicators[-match_index,] <- gamma_nonmatch
  }
  Sadinle_indicators <- purrr::map2(data.frame(indicators), levels, ~fs_to_sadinle_2(.x, .y)) %>%
    do.call(cbind, .)

  list(comparisons = Sadinle_indicators,
       n1 = n1,
       n2 = n2,
       nDisagLevs = levels,
       Ztrue = Ztrue)
}

simulate_comparisons_mm <- function(m, u, levels, n1, n2, overlap){
  field_marker <- unlist(lapply(1:length(levels), function(x){
    rep(x, levels[x])
  }))

  N <- n1 * n2
  ids <- expand.grid(1:n1, 1:n2)
  indicators <- matrix(NA, nrow = N, ncol = length(levels))



  df2_matches <- lapply(overlap, seq_len)
  matches_in_2 <- seq_len(overlap[1])
  match_breaks <- cumsum(c(0, overlap[ -length(overlap)]))
  df1_matches <- purrr::map2(df2_matches, match_breaks, ~.x + .y)

  Z_temp <- data.frame(match_1 = unlist(df1_matches),
                       match_2 = unlist(df2_matches))

  # Z_true <- data.frame(match_2 = 1:n2, match_1 = 0) %>%
  #   nest_by(match_2, .key = "match_1")
  #
  # Z_matches <- Z_temp %>%
  #   nest_by(match_2, .key = "match_1", .keep = F)
  #
  # Z_true[matches_in_2, ]$match_1 <- Z_matches$match_1

  match_index <- apply(Z_temp, 1, function(x){
    which(ids[, 1] == unlist(x)[1] & ids[, 2] == unlist(x)[2])
  })

  Z_true <- Z_temp %>%
    setNames(c("id_1", "id_2")) %>%
    tidyr::complete(id_2 = 1:n2) %>%
    relocate(id_2, .after = id_1)


  # %>%
  #   nest_by(id_2)


  m_list <- split(m, field_marker)
  u_list <- split(u, field_marker)

  gamma_match <- sapply(m_list, function(x){
    sample(seq_along(x) - 1, sum(overlap), replace = T, x)
  })

  gamma_nonmatch <- sapply(u_list, function(x){
    sample(seq_along(x) - 1, N - sum(overlap), replace = T, x)
  })

  if(sum(overlap) == 0){
    indicators <- gamma_nonmatch
  } else {
    indicators[match_index,] <- gamma_match
    indicators[-match_index,] <- gamma_nonmatch
  }
  Sadinle_indicators <- purrr::map2(data.frame(indicators), levels, ~fs_to_sadinle_2(.x, .y)) %>%
    do.call(cbind, .)


  list(comparisons = Sadinle_indicators,
       n1 = n1,
       n2 = n2,
       nDisagLevs = levels,
       Z_true = Z_true)
}

