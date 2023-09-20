hash_comparisons__mm<- function(cd,
                    method = "both", R = NULL,
                    all_patterns = FALSE){


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

  temp2 <- temp %>%
    group_split(rec2)

  thing <- temp2[[1]]

  thing$hash_id %>%
    as.numeric() %>%
    combn(2) %>%
    t() %>%
    data.frame() %>%
    group_by(X1, X2) %>%
    count(n)

  thing$hash_id %>%
    as.numeric() %>%
    combn(1) %>%
    t() %>%
    data.frame() %>%
    group_by_all() %>%
    count()

  lookat <- combn(as.numeric(thing$hash_id), 2, simplify = F)
