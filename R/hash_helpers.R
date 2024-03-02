#'
#' Internal Functions
#'
#'
#'
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
