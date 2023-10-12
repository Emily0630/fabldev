#' @export
#'
combine_hash <- function(hash_list, n1, n2){

    total_counts <- purrr::map(hash_list, ~.x$total_counts) %>%
    do.call(rbind, .) %>%
    colSums()

  pattern_counts_by_record <- hash_list %>%
    purrr::map(`[[`, "pattern_counts_by_record") %>%
    purrr::flatten()

  hash_to_file_1 <- hash_list %>%
    purrr::map(`[[`, "hash_to_file_1") %>%
    purrr::flatten()

  record_counts_by_pattern <- purrr::transpose(pattern_counts_by_record) %>%
    purrr::map(unlist) %>%
    purrr::map(unname)

  flags <- hash_list %>%
    purrr::map(`[[`, "flags") %>%
    purrr::flatten()

  patterns <- list(ohe = hash_list[[1]]$ohe,
                   total_counts = total_counts,
                   pattern_counts_by_record = pattern_counts_by_record,
                   record_counts_by_pattern = record_counts_by_pattern,
                   hash_to_file_1 = hash_to_file_1,
                   flags = flags,
                   field_marker = hash_list[[1]]$field_marker,
                   n1 = n1,
                   n2 = n2)

}
