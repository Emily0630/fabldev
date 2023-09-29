hash_comparisons__mm<- function(cd,
                    method = "both", R = NULL,
                    all_patterns = TRUE,
                    max_K = 1){


  indicators <- cd[[1]]
  N <- dim(indicators)[1]
  n1 <- cd[[2]]
  n2 <- cd[[3]]

  levels <- cd[[4]]
  fields <- seq_along(cd[[4]])
  field_marker <- sapply(fields, function(x){
    rep(x, cd[[4]][x])
  }) %>%
    unlist(.) %>%
    as.vector(.)

  ids <- expand.grid(1:n1, 1:n2)
  rec1 <- ids[,1]
  rec2 <- ids[,2]

  Lf_vec<- (levels) %>%
    c(0, .) %>%
    cumsum()

  hash_vals <- purrr::imap(cd[[4]], ~hash_field(.x, .y, Lf_vec)) %>%
    unlist()

  hash <- sweep(indicators, 2, hash_vals, "*") %>%
    rowSums() + 1

  if(all_patterns == TRUE){

    unique_patterns <- possible_patterns_sadinle(levels)
    unique_hashed <- sweep(unique_patterns, 2, hash_vals, "*") %>%
      rowSums() + 1
    P <- dim(unique_patterns)[1]
    hash_id <- match(hash, unique_hashed) %>%
      factor(levels = 1:P)

  } else {

    unique_hashed <- unique(hash)
    P <- length(unique_hashed)
    hash_id <- match(hash, unique_hashed) %>%
      factor(levels = 1:P)
    unique_patterns <- indicators[!duplicated(hash_id), ]
  }


  temp <- data.frame(indicators, rec1, rec2, hash_id)

  max_K <- 2
  column_names <- sapply(seq_len(max_K), function(x){
    paste0("X", x)
  })
  # pattern_combinations <- lapply(seq_len(max_K), function(x){
  #   if(x == 1){
  #     seq(1, P)
  #   } else {
  #     seq(0, P)
  #   }
  # }) %>%
  #   do.call(expand.grid, .) %>%
  #   data.frame() %>%
  #   setNames(column_names)

  pattern_combinations <- lapply(seq_len(max_K), function(k){
    do.call(expand.grid, rep(list(1:P), k)) %>%
      data.frame() %>%
      setNames(column_names[1:k])
  })

  combination_ohe <- lapply(seq_len(max_K), function(k){
    apply(pattern_combinations[[k]], 1, function(x){
      if(k == 1){
        unique_patterns[x, ]
      } else {
      unique_patterns[unlist(x), ] %>%
        colSums()
      }
    }) %>%
      t()
  })

  temp2 <- temp %>%
    group_split(rec2)

  combination_counts <- lapply(temp2, function(y){
    combo_counts <- lapply(seq_len(max_K), function(k){
      combo_counts <- y$hash_id %>%
        as.numeric() %>%
        combn(k) %>%
        t() %>%
        data.frame() %>%
        setNames(column_names[1:k]) %>%
        group_by_all() %>%
        count()

      full_set <- left_join(pattern_combinations[[k]], combo_counts)
      full_set$n[is.na(full_set$n)] <- 0
      full_set$n
    })
    combo_counts
  })

  hash_to_file_1 <- temp %>%
    select(rec1, rec2, hash_id) %>%
    nest_by(rec2, hash_id, .keep = F) %>%
    mutate(hash_id = as.integer(hash_id)) %>%
    rowwise() %>%
    mutate(N = nrow(data))

  hash_to_file_1 <- left_join(x = pattern_lookup,
                              y = hash_to_file_1,
                              by = c("hash_id", "rec2"))



