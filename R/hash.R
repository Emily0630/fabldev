hash_field <- function(L_f, k, Lf_vec){
  level_seq <- seq_len(L_f)
  as.numeric(level_seq > 0) * 2 ^ ((level_seq) + (as.numeric(k > 1)  * Lf_vec[k]))
}

hash_comparisons <- function(comparisons,
                 R = NULL,
                 all_patterns = FALSE){

  indicators <- comparisons[[1]]
  n1 <- comparisons[[2]]
  n2 <- comparisons[[3]]

  levels <- comparisons[[4]]
  fields <- seq_along(comparisons[[4]])
  field_marker <- sapply(fields, function(x){
    rep(x, comparisons[[4]][x])
  }) %>%
    unlist(.) %>%
    as.vector(.)

  ids <- expand.grid(1:n1, 1:n2)
  rec1 <- ids[,1]
  rec2 <- ids[,2]
  levels <- cd[[4]]

  Lf_vec<- (levels) %>%
    c(0, .) %>%
    cumsum()

  hash_vals <- purrr::imap(cd[[4]], ~hash_field(.x, .y, Lf_vec)) %>%
    unlist()

  hash <- sweep(indicators, 2, hash_vals, "*") %>%
    rowSums() + 1

  if(all_patterns == TRUE){
    unique_patterns <- GetPossiblePatternsSad_missing(levels)
    unique_hashed <- sweep(unique_patterns, 2, hash_vals, "*") %>%
      rowSums() + 1
    P <- dim(unique_patterns)[1]
    hash_id <- match(hash, unique_hashed) %>%
      factor(., 1:P)

  } else {

    unique_hashed <- unique(hash)
    P <- length(unique_hashed)
    hash_id <- match(hash, unique_hashed) %>%
      factor(., 1:P)
    unique_patterns <- indicators[!duplicated(hash_id),]
  }


  temp <- data.frame(indicators, rec1, rec2, hash_id)
  pattern_counts <- temp %>%
    group_by(hash_id, .drop = F) %>%
    count() %>%
    pull()

  thing <- split(temp, rec2)
  hash_to_file_1 <- purrr::map(thing, ~ lapply(1:P, function(y){
    which(.x$hash_id == y)
  })
  )

  counts_by_rec <-  purrr::map(hash_to_file_1, ~lapply(.x, length)) %>%
    purrr::map(unlist)

  hash_to_file_1 <- lapply(hash_to_file_1, function(x){
    append(x, list(n1 +1))
  })

  if(!is.null(R)){
    hash_to_file_1 <- lapply(hash_to_file_1, function(z){
      purrr::map(z, ~sei(.x, R))
    })
  }


  patterns <- list(ohe = unique_patterns,
                   total_counts = pattern_counts,
                   pattern_counts_by_record = counts_by_rec,
                   hash_to_file_1 = hash_to_file_1,
                   field_marker = field_marker,
                   n1 = n1,
                   n2 = n2)
  patterns

}
