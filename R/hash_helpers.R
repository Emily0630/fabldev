#'
#' Internal Functions
#'
#'
#'
hash_field <- function(L_f, k, Lf_vec){
  level_seq <- seq_len(L_f)
  as.numeric(level_seq > 0) * 2 ^ ((level_seq) + (as.numeric(k > 1)  * Lf_vec[k]))
}

fs_to_ohe <- function(gamma, levels){
  unname(unlist(mapply(function(z, y){
    gamma_f <- rep(0, y)
    if(z == 0){
      return(gamma_f)
    }
    gamma_f[z] <- 1
    gamma_f
  }, z = gamma, y = levels)))
}

possible_patterns_ohe <- function(levels){
  possible_values <- lapply(levels, function(x){
    c(0, seq_len(x))
  })
  possible_patterns <- data.frame(do.call(expand.grid, possible_values))

  data.frame(t(apply(possible_patterns, 1, function(x){
    fs_to_ohe(x, levels)
  })))

}

sei <- function(x, R){
  if(length(x) <= R){
    return(x)
  } else {
    sample(x, R, replace = FALSE)
  }
}
