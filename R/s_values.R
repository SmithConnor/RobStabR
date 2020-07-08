#' Compute the subtractive lack-of-fit measures
#'
#' @param weights
#' @param data a data frame containing the variables in the model.
#' @param n
#' @param p
#' @param family
#' @param coef
#' @param wald
#' @param dev
#' @importFrom magrittr %>%
#' @export
#'
#' @example

s_values = function(weights,
                    data,
                    n,
                    p,
                    family,
                    coef,
                    wald,
                    dev,
                    tcc){
  glmBoot = robustbase::glmrob(formula = y~.,
                               data = data,
                               subset = weights,
                               family = family,
                               control = robustbase::glmrobMqle.control(tcc = tcc, maxit = 1000))
  glmBootSummary = base::summary(glmBoot)
  sVal = base::matrix(data = NA_integer_,
                      nrow = 3,
                      ncol = p)
  QD = base::rep(x = NA,
                 times = p)
  colnames(sVal) = base::names(glmBoot$coefficients[-1])
  rownames(sVal) = c("Coef", "Wald", "Dev")
  if(coef == TRUE){
    sVal[1,] = glmBootSummary$coefficients[-1, 1]
  }
  if(wald == TRUE){
    sVal[2,] = glmBootSummary$coefficients[-1, 3]
  }
  if(dev == TRUE){
    for(i in 1:p){
      variables = utils::head(x = names(data),
                              n = p)
      formulaDev = as.formula(paste( "y ~", paste(variables[-i], collapse = "+")))
      glmDev = robustbase::glmrob(formula = formulaDev,
                                  data = data,
                                  subset = weights,
                                  family = family,
                                  control = robustbase::glmrobMqle.control(tcc = tcc, maxit = 1000))
      anovaDev = base::tryCatch(list(fit = anova(glmDev,
                                                 glmBoot,
                                                 test = "QD"),
                                     QD = TRUE),
                                error = function(e){
                                  'NA'
                                }

      )

      if(base::all(anovaDev == 'NA')){
        anovaDev = list(fit = anova(glmDev,
                                    glmBoot,
                                    test = "QDapprox"),
                        QD = FALSE)
      }

      sVal[3,i] = anovaDev$fit$Test.Stat[2]
      QD[i] = anovaDev$QD
    }
  }
  return(list(sVal, QD))
}

#####

#' Compute the subtractive lack-of-fit measures
#'
#' @param list
#' @param row
#'
#' @example

solution_path = function(list,
                         row){
  matrix = list[[1]]
  vector = matrix[row, ]
  p = base::length(vector)
  absVector = base::abs(vector)
  ranks =  p + 1 - base::rank(absVector)
  var = base::names(ranks)
  n = base::length(ranks)
  matrix = base::matrix(data = NA_integer_,
                        nrow = n-1,
                        ncol = 2)
  df = base::data.frame(matrix)
  base::names(df) <- c("Variables","Dimension")
  for (i in 1:(n-1)){
    whichVar = ranks <= i
    string = base::paste(c("1", var[whichVar]), sep = "", collapse = "+")
    df$Variables[i] = string
    df$Dimension[i] = i
  }
  return(df)
}
