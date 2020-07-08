#'
#' @param data
#' @param B
#' @param m
#' @param nStrata
#' @param family
#' @param k
#' @param resid
#' @param coef
#' @param wald
#' @param dev
#' @export
#'
#' @example



model_space_non = function(data,
                       B,
                       m,
                       nStrata = 8,
                       family,
                       k = 1,
                       resid = "pearson",
                       coef = TRUE,
                       wald = TRUE,
                       dev = TRUE,
                       bootstraps = NA){
  # Parameters
  n = base::NROW(data)
  p = base::NCOL(data) - 1
  nMethods = coef + wald + dev
  # Variables
  bootstrapModels = base::rep(x = list(NA),
                              times = B)
  varNames = base::colnames(x = data)[-(p+1)]
  coefValues = base::matrix(data = NA_integer_,
                            nrow = B,
                            ncol = p)
  colnames(coefValues) = varNames
  waldValues = matrix(data = NA_integer_,
                      nrow = B,
                      ncol = p)
  colnames(waldValues) = varNames
  devValues = matrix(data = NA_integer_,
                     nrow = B,
                     ncol = p)
  colnames(devValues) = varNames
  output = base::list()
  sVals = base::list()
  # Initial Fit
  glmFull = stats::glm(formula = y~.,
                       data = data,
                       family = family)
  glmFullResid = stats::residuals(object = glmFull,
                                  type = resid)
  residOrder = rank(x = glmFullResid)
  residStrata = ceiling(x = residOrder * nStrata / n)
  # Generate weights
  if(base::all(base::is.na(bootstraps))){
    bootstraps = base::matrix(data = NA_integer_,
                              nrow = B,
                              ncol = m)
    for(b in 1:B){
      for(o in 1:nStrata){
        bootstraps[b, (o - 1) * m / nStrata + (1 : (m / nStrata))] = base::sample(x = which(residStrata == o),
                                                                                  size = m / nStrata)
      }
    }}
  # Fit Models
  sValues = plyr::alply(.data = bootstraps,
                        .margins = 1,
                        .fun = s_values_non,
                        data = data,
                        n= n,
                        p = p,
                        family = family,
                        coef = coef,
                        wald = wald,
                        dev = dev,
                        .parallel = FALSE)
  if(coef == TRUE){
    coefPath = lapply(X = sValues,
                      FUN = solution_path_non,
                      row = 1)
    coefSpace = reduced_space(list = coefPath, vector = varNames)
    output$coefSpace = coefSpace
    sVals$coef = base::matrix(data = NA_integer_,
                              nrow = B,
                              ncol = p)
    for (i in 1:B){
      sVals$coef[i,] = sValues[[i]][1,]
    }
    base::colnames(sVals$coef) = varNames
    base::rownames(sVals$coef) = base::paste0("Bootstrap",1:B)
  }
  if(wald == TRUE){
    waldPath = lapply(X = sValues,
                      FUN = solution_path_non,
                      row = 2)
    waldSpace = reduced_space(list = waldPath, vector = varNames)
    output$waldSpace = waldSpace
    sVals$wald = base::matrix(data = NA_integer_,
                              nrow = B,
                              ncol = p)
    for (i in 1:B){
      sVals$wald[i,] = sValues[[i]][2,]
    }
    base::colnames(sVals$wald) = varNames
    base::rownames(sVals$wald) = base::paste0("Bootstrap",1:B)
  }
  if(dev == TRUE){
    devPath = lapply(X = sValues,
                     FUN = solution_path_non,
                     row = 3)
    devSpace = reduced_space(list = devPath, vector = varNames)
    output$devSpace = devSpace
    sVals$dev = base::matrix(data = NA_integer_,
                             nrow = B,
                             ncol = p)
    for (i in 1:B){
      sVals$dev[i,] = sValues[[i]][3,]
    }
    base::colnames(sVals$dev) = varNames
    base::rownames(sVals$dev) = base::paste0("Bootstrap",1:B)
  }

  countModels = count_models(output,
                             k = k)
  additional = countModels %>%
    base::names()
  addition = base::paste0(additional, "Count")
  for(i in 1:nMethods){
    output[[nMethods + i]] = countModels[[i]]
  }
  names(output)[nMethods + 1:nMethods] = addition
  for(i in 1:nMethods){
    output[[2*nMethods + i]] = sVals[[i]]
  }
  names(output)[2*nMethods + 1:(nMethods)] = base::paste0(additional, "SVal")
  output$bootsrap = bootstraps
  output$varNames = varNames
  return(output)
}


###############

#'
#' @param weights
#' @param data
#' @param n
#' @param p
#' @param family
#' @param coef
#' @param wald
#' @param dev
#'
#' @example

s_values_non = function(weights,
                    data,
                    n,
                    p,
                    family,
                    coef,
                    wald,
                    dev){
  glmBoot = stats::glm(formula = y~.,
                       data = data,
                       subset = weights,
                       family = family)
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
      glmDev = stats::glm(formula = formulaDev,
                          data = data,
                          subset = weights,
                          family = family)
      anovaDev = anova(glmDev, glmBoot, test = "Chisq")

      sVal[3,i] = anovaDev$Deviance[2]
    }
  }
  return(sVal)
}


####################

IC_non = function(vector, data, family){
  glmFull = stats::glm(y~.,
                       data = data,
                       family = family)
  model = base::paste0("y~", vector[1]) %>%
    stats::as.formula(.)
  glmFit = stas::glm(formula = model,
                     data = data,
                     family = family)
  dev = anova(glmFit, glmFull, test = "Chisq")$Deviance[2]
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

solution_path_non = function(list, row){
  matrix = list
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
