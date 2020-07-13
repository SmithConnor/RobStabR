#' Calculate the Information Criterion for a robust model
#'
#' @param vector a boolean vector indicating which variables to include.
#' @param data a data frame containing the variables in the model.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param tcc the tuning constant c in Huber's psi-function.
#' @importFrom magrittr %>%

IC = function(vector,
              data,
              family,
              tcc){
  glmFull = robustbase::glmrob(y~.,
                               data = data,
                               family = family,
                               control = robustbase::glmrobMqle.control(tcc = tcc, maxit = 1000))
  model = base::paste0("y~", vector[1]) %>%
    stats::as.formula()
  glmFit = robustbase::glmrob(formula = model,
                              data = data,
                              family = family,
                              control = robustbase::glmrobMqle.control(tcc = tcc, maxit = 1000))
  dev = stats::anova(glmFit, glmFull, test = "QDapprox")$Test.Stat[2]
  ###
  n = nrow(data)
  p = base::length(glmFit$coefficients) - 1
  ###
  out = base::data.frame(model = vector[1],
                         dimension = vector[2],
                         count = vector[3],
                         P1 = as.numeric(dev  + p*log(n)),
                         P2 = as.numeric(dev + p*(log(n) + 1)))
  return(out)
}


#####

#' Calculate the Information Criterion for a robust model over a matrix
#'
#' @param matrix a boolean matrix indicating which variables to include on each row.
#' @param data a data frame containing the variables in the model.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param tcc the tuning constant c in Huber's psi-function.
#' @importFrom magrittr %>%

matrix_IC = function(matrix,
                     data,
                     family,
                     tcc){
  tictoc::tic()
  output = base::apply(X = matrix,
              MARGIN = 1,
              FUN = IC,
              data = data,
              family = family,
              tcc = tcc)
  endTime = tictoc::toc(quiet = TRUE)
  return(base::list(candidateSpace = output,
         timing = endTime$toc - endTime$tic))

}
