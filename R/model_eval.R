#' @param df
#' @param criteria
#' @importFrom magrittr %>%

best_model = function(df,
                      criteria){
  lowest = df[,criteria] %>%
    which.min()
  df$Best = FALSE
  df$Best[lowest] = TRUE
  return(df)
}

#####

#' @param RSSM
#' @param penalty
#' @importFrom magrittr %>%

best_robust = function(RSSM,
                       penalty){
  coefBest = RSSM$coef[,penalty] %>% base::which.min(.)
  coefOut = RSSM$coef$model[coefBest] %>%
    as.character() %>%
    base::strsplit(., "\\+")
  coefOut = coefOut[[1]][-1]
  waldBest = RSSM$wald[,penalty] %>% base::which.min(.)
  waldOut = RSSM$wald$model[waldBest] %>% as.character() %>%
    base::strsplit(., "\\+")
  waldOut = waldOut[[1]][-1]
  devBest = RSSM$dev[,penalty] %>% base::which.min(.)
  devOut = RSSM$dev$model[devBest] %>% as.character() %>%
    base::strsplit(., "\\+")
  devOut = devOut[[1]][-1]
  return(list(coef = coefOut, wald = waldOut, dev = devOut))
}

#####

#'
#' @param data a data frame containing the variables in the model.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param anovaTest
#' @param pVal
#' @importFrom magrittr %>%
#' @import formula.tools

step_glmrob = function(data,
                       family,
                       anovaTest = "QDapprox",
                       pVal = 0.05,
                       tcc){
  glmFull = robustbase::glmrob(y~.,
                               data = data,
                               family = family,
                               control = robustbase::glmrobMqle.control(tcc = tcc, maxit = 1000))
  glmSummary = list()
  variables = names(data)
  p = ncol(data) - 1
  variables = head(x = variables,
                   n = p)
  vars = variables
  glmSteps = matrix(data = NA,
                    nrow = 3,
                    ncol = p)
  base::row.names(glmSteps) = c("Initial Model", "Remove variable", "p-value")
  continue = TRUE
  iter = 1
  mCount = 1
  while(continue == TRUE){
    pval = base::rep(x = 0,
                     times = length(vars))
    model = base::rep(x = NA,
                      times = length(vars))
    for(k in 1:length(vars)){
      formula = stats::as.formula(base::paste( "y ~", paste(vars[-k], collapse = "+")))
      glmFit = robustbase::glmrob(formula,
                                  data = data,
                                  family = family,
                                  control = robustbase::glmrobMqle.control(tcc = tcc, maxit = 1000))
      mCount = mCount + 1
      anovaFit = anova(glmFit,
                       glmFull,
                       test = "QDapprox")
      model[k] = as.character(x = formula)
      pval[k] = anovaFit$`Pr(>chisq)`[2]
    }
    check = sum(pval >= pVal)
    if(check == 0){
      continue = FALSE
    }
    remove = base::which.max(x = pval)
    glmSteps[,iter] = c(base::paste( "y ~", base::paste(vars, collapse = "+")), vars[remove] , pval[remove])
    vars = vars[-remove]
    iter = iter + 1
  }
  glmSteps = janitor::remove_empty(glmSteps, which = "cols")
  base::colnames(glmSteps) = base::paste0("Step", 1:base::NCOL(glmSteps))
  glmSteps = glmSteps %>%
    data.frame(., stringsAsFactors = FALSE)
  final = base::substring(glmSteps[1,(iter - 1)], 5) %>%
    base::strsplit(., "\\+")
  final = final[[1]]
  glmSteps = glmSteps[,-c(iter - 1)]
  return(final)
}

#####

#'
#' @param data a data frame containing the variables in the model.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param k
#' @importFrom magrittr %>%

step_ic = function(data,
                   family,
                   k = 2){
  glmFull = stats::glm(y~.,
                       data = data,
                       family = family)
  step = MASS::stepAIC(object = glmFull,
                       direction = "backward",
                       trace = 0,
                       k = k)
  stepModel = step$model %>%
    base::colnames()
  stepModel = stepModel[-1]
  return(stepModel)
}
