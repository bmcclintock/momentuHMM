# check that the formula is a formula
is.formula <- function(x)
  tryCatch(inherits(x,"formula"),error= function(e) {FALSE})