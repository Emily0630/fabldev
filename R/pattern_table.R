#' @export
#'
fs_to_sadinle <- function(gamma, levels){
  unname(unlist(mapply(function(z, y){
    gamma_f <- rep(0, y)
    if(z == 0){
      return(gamma_f)
    }
    gamma_f[z] <- 1
    gamma_f
  }, z = gamma, y = levels)))
}

fs_to_sadinle_2 <- function(gamma_f, Lf){
  new_gamma <- matrix(0, nrow = length(gamma_f), ncol = Lf)
  for(i in seq_along(gamma_f)){
    new_gamma[i, gamma_f[i]+1] <- 1
  }
  new_gamma
}

possible_patterns_sadinle <- function(levels){
  possible_values <- lapply(levels, function(x){
    c(0, seq_len(x))
  })
  possible_patterns <- data.frame(do.call(expand.grid, possible_values))

  thing <- data.frame(t(apply(possible_patterns, 1, function(x){
    fs_to_sadinle(x, levels)
  })))
  thing
}
