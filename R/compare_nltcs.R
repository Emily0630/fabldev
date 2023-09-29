
compare_binary <- function(df1, df2){
  nA <- length(df1)
  nB <- length(df2)[1]
  ids <- expand.grid(1:nA, 1:nB)
  match <- df1[ids[, 1]]== df2[ids[, 2]]
  comparison <- cbind(match, 1 - match) %>%
    unname()
  comparison
}

compare_dob <- function(vec1, vec2){
  nA <- length(vec1)
  nB <- length(vec2)
  ids <- expand.grid(1:nA, 1:nB)
  matched_components <- matrix(NA, nrow = nrow(ids), ncol = 3)
  for(i in 1:3){
    s <- 2*i - 1
    match <- stringr::str_sub(vec1[ids[, 1]], s, s+1) == stringr::str_sub(vec2[ids[, 2]], s, s+1)
    matched_components[, i] <- match
  }
  score <- rowSums(matched_components)
  score[score == 0] <- 1

  comparison <- matrix(0, nrow = nrow(ids), ncol = 3)
  comparison[score == 3, 1] <- 1
  comparison[score == 2, 2] <- 1
  comparison[score == 1, 3] <- 1
  comparison
}

compare_location <- function(vec1, vec2){
  nA <- length(vec1)
  nB <- length(vec2)
  ids <- expand.grid(1:nA, 1:nB)
  matched_components <- matrix(NA, nrow = nrow(ids), ncol = 2)
  for(i in 1:2){
    s <- 2*i - 1
    match <- stringr::str_sub(vec1[ids[, 1]], s, s+1) == stringr::str_sub(vec2[ids[, 2]], s, s+1)
    matched_components[, i] <- match
  }
  score <- rep(1, nrow(ids))
  score[matched_components[, 2] == 1] <- 3
  score[matched_components[, 2] == 1 & matched_components[, 1] == 0] <- 2

  comparison <- matrix(0, nrow = nrow(ids), ncol = 3)
  comparison[score == 3, 1] <- 1
  comparison[score == 2, 2] <- 1
  comparison[score == 1, 3] <- 1
  comparison
}

# compare_nltcs <- function(df1, df2){
#   nA <- nrow(df1)
#   nB <- nrow(df2)
#   cd1 <- compare_binary(df1$SEX, df2$SEX)
#   cd2 <- compare_dob(df1$dob, df2$dob)
#   cd3 <- compare_location(df1$location, df2$location)
#   comparison <- cbind(cd1, cd2, cd3)
#
#   cd <- list(comparisons = comparison,
#              n1 = nA,
#              n2 = nB,
#              nDisagLevs = c(2, 3, 3),
#              Ztrue = NULL)
#   cd
# }

compare_nltcs <- function(df1, df2){
  nA <- nrow(df1)
  nB <- nrow(df2)
  cd1 <- compare_binary(df1$SEX, df2$SEX)
  cd2 <- compare_binary(df1$DOB_DAY, df2$DOB_DAY)
  cd3 <- compare_binary(df1$DOB_MONTH, df2$DOB_MONTH)
  cd4 <- compare_binary(df1$DOB_YEAR, df2$DOB_YEAR)
  cd5 <- compare_location(df1$location, df2$location)
  comparison <- cbind(cd1, cd2, cd3, cd4, cd5)

  cd <- list(comparisons = comparison,
             n1 = nA,
             n2 = nB,
             nDisagLevs = c(2, 2, 2, 2, 3),
             Ztrue = NULL)
  cd
}
