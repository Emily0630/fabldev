hash_comparisons_mm<- function(cd,
                    method = "both", R = NULL,
                    all_patterns = FALSE,
                    max_K){


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
  total_counts <- temp %>%
    group_by(hash_id, .drop = F) %>%
    count() %>%
    pull()

  column_names <- sapply(seq_len(max_K), function(x){
    paste0("X", x)
  })

  pattern_lookup <- expand.grid(1:P, 1:n2) %>%
    data.frame() %>%
    setNames(., c("hash_id", "rec2"))

  pattern_combinations <- lapply(seq_len(max_K), function(k){
    do.call(expand.grid, rep(list(1:P), k)) %>%
      data.frame() %>%
      setNames(column_names[1:k])
  })

  pattern_combinations_all <- do.call(plyr::rbind.fill, pattern_combinations) %>%
    data.frame() %>%
    mutate(rn = row_number())

  K_marker <- apply(pattern_combinations_all, 1, function(x){
    sum(!is.na(x)) - 1
  })

  pattern_finder <- lapply(1:P, function(y){
    pattern_combinations_all %>%
      filter(if_any(starts_with("X"), ~.x == y )) %>%
      select(rn) %>%
      pull()
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

  # by_names <- sapply(1:max_K, function(x){
  #   paste0("X", x)
  # })

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

      by_names <- sapply(1:k, function(x){
        paste0("X", x)
      })

      full_set <- left_join(pattern_combinations[[k]],
                            combo_counts,
                            by = by_names)
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

  hash_to_file_1$N[is.na(hash_to_file_1$N)] <- 0

  # pattern_counts_by_record <- hash_to_file_1 %>%
  #   select(-data) %>%
  #   group_split(rec2, .keep = F) %>%
  #   purrr::map(., `[[`, "N")

  pattern_counts_by_record <- lapply(seq_len(max_K), function(k){
    combination_counts %>%
      purrr::map(., `[[`, k)
  })


  # pattern_counts_by_record <- combination_counts %>%
  #   purrr::map(., `[[`, 1)

  record_counts_by_pattern <- lapply(seq_len(max_K), function(k){
    purrr::transpose(pattern_counts_by_record[[k]]) %>%
      purrr::map(unlist) %>%
      purrr::map(unname)
  })

  # total_counts_K <- lapply(seq_len(max_K), function(k){
  # record_counts_by_pattern[[k]] %>%
  #   do.call(cbind, .) %>%
  #   colSums()
  # })

  # record_counts_by_pattern <- purrr::transpose(pattern_counts_by_record) %>%
  #   purrr::map(unlist) %>%
  #   purrr::map(unname)

  hash_to_file_1 <- hash_to_file_1 %>%
    group_split(rec2) %>%
    purrr::map(., ~ .x %>%
                 group_split(hash_id)) %>%
    purrr::map(., ~purrr::map(.x, `[[`, "data")) %>%
    purrr::map(., ~purrr::map(., ~ unname(unlist(.x))))

  if(!is.null(R)){
    hash_to_file_1 <- lapply(hash_to_file_1, function(z){
      purrr::map(z, ~sei(.x, R))
    })}


  patterns <- list(ohe = unique_patterns,
                   combination_ohe = combination_ohe,
                   total_counts = total_counts,
                   combination_counts = combination_counts,
                   pattern_counts_by_record = pattern_counts_by_record,
                   record_counts_by_pattern = record_counts_by_pattern,
                   hash_to_file_1 = hash_to_file_1,
                   #flags = flags,
                   field_marker = field_marker,
                   n1 = n1,
                   n2 = n2,
                   max_K = max_K,
                   #total_counts_K = total_counts_K,
                   pattern_finder = pattern_finder,
                   K_marker = K_marker)
}



