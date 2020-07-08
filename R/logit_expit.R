#' Expit function
#' @param x numeric
#' @export
#' @examples
#' curve(expit,
#'  from = -10,
#'  to = 10)

expit = function(x){
  return(1/(1+exp(-x)))
}

#####

#' Logit function
#' @param x numeric
#' @export
#' @examples
#' curve(logit,
#' from = 0.01,
#' to = 0.99)

logit = function(x){
  if(any(x <= 0) | any(x >= 1)){
    error("x must be between zero and 1")
  }
  return(log(x) - log(1-x))
}
