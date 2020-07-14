#' Variable selection using the bestglm package
#'
#' @param data a data frame containing the variables in the model.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param IC information criteria to use: "AIC", "BIC", "BICg", "BICq", "LOOCV", "CV".

bestglm_ic = function(data,
                      family,
                      IC = "AIC"){
  search = bestglm::bestglm(Xy = data,
                            family = family,
                            IC = IC)
  searchOut = names(search$BestModel$coefficients)[-1]
  return(searchOut)
}

#####

#' Exhaustive variable selection using the RDBC
#'
#' @param data a data frame containing the variables in the model.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param tcc the tuning constant c in Huber's psi-function.

exhaustive_RDBC = function(data,
                           family,
                           tcc){
  n = base::nrow(data)
  p = base::ncol(data) - 1
  l = base::rep(list(0:1), p)
  expl = expand.grid(l)
  varNames = base::colnames(data)[-(p+1)]
  combForm = function(vector, names){
    stats::as.formula(base::paste0("y~",
                                   base::paste0(c("1",varNames[base::as.logical(vector)]), collapse = "+"),
                                   collapse=""))
  }
  allForm = base::apply(expl, 1, combForm, names = varNames)
  glmFit = lapply(allForm,
                  robustbase::glmrob,
                  data = data,
                  family = family,
                  control = robustbase::glmrobMqle.control(tcc = tcc, maxit = 1000))
  anovaFit = list()
  for(i in 1:(length(glmFit)-1)){
    anovaFit[[i]] = anova(glmFit[[i]], glmFit[[2^p]], test = "QDapprox")$Test.Stat[2]
  }
  RDCB = base::unlist(anovaFit) + apply(expl,1,sum)[-(2^p)]*(log(n) + 1)
  final = varNames[as.logical(expl[base::which.min(RDCB),])]
  return(final)
}
