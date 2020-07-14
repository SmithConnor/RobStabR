#' Find the candidate model space after B bootstraps
#'
#' @param data a data frame containing the variables in the model.
#' @param B The number of bootstrap re-samples to run.
#' @param m The m used for the m-out-of-n boostrap.
#' @param nStrata Number of strat used for the stratified re-sampling.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param k The minimum number of times a mondel must be identified through resampling to be included in the rediced candidate space.
#' @param resid the type of residuals to be used for the stratified re-sampling.
#' @param coef A TRUE/FALSE value to indicate whether to evaluate RobStab using regression coefficients.
#' @param wald A TRUE/FALSE value to indicate whether to evaluate RobStab using wald statisitcs.
#' @param dev A TRUE/FALSE value to indicate whether to evaluate RobStab using deviances.
#' @param bootstraps Recreact results using previously used resamples.
#' @param tcc the tuning constant c in Huber's psi-function.
#' @importFrom magrittr %>%
#' @export
#'



model_space = function(data,
                       B,
                       m,
                       nStrata = 8,
                       family,
                       k = 1,
                       resid = "pearson",
                       coef = TRUE,
                       wald = TRUE,
                       dev = TRUE,
                       bootstraps = NA,
                       tcc){
  tictoc::tic()
  # Parameters
  n = base::NROW(data)
  p = base::NCOL(data) - 1
  nMethods = coef + wald + dev
  # Variables
  devQD = base::matrix(0, nrow = B, ncol = p)
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
  glmFull = robustbase::glmrob(formula = y~.,
                               data = data,
                               family = family,
                               control = robustbase::glmrobMqle.control(tcc = tcc,
                                                                        maxit = 1000))
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
  endStart = tictoc::toc(quiet = TRUE)
  timing = endStart$toc - endStart$tic
  # Fit Models
  sValues = plyr::alply(.data = bootstraps,
                        .margins = 1,
                        .fun = s_values,
                        data = data,
                        n= n,
                        p = p,
                        family = family,
                        coef = coef,
                        wald = wald,
                        dev = dev,
                        .parallel = FALSE,
                        tcc = tcc)
  timing3 = base::rep(0,
                      times = 3)
  if(coef == TRUE){
    tictoc::tic()
    coefPath = lapply(X = sValues,
                      FUN = solution_path,
                      row = 1)
    coefSpace = reduced_space(list = coefPath, vector = varNames)
    output$coefSpace = coefSpace
    sVals$coef = base::matrix(data = NA_integer_,
                              nrow = B,
                              ncol = p)
    for (i in 1:B){
      sVals$coef[i,] = sValues[[i]][[1]][1,]
    }
    base::colnames(sVals$coef) = varNames
    base::rownames(sVals$coef) = base::paste0("Bootstrap",1:B)
    coefTime = tictoc::toc(quiet = TRUE)
    timing3[1] = coefTime$toc - coefTime$tic
  }
  if(wald == TRUE){
    tictoc::tic()
    waldPath = lapply(X = sValues,
                      FUN = solution_path,
                      row = 2)
    waldSpace = reduced_space(list = waldPath, vector = varNames)
    output$waldSpace = waldSpace
    sVals$wald = base::matrix(data = NA_integer_,
                              nrow = B,
                              ncol = p)
    for (i in 1:B){
      sVals$wald[i,] = sValues[[i]][[1]][2,]
    }
    base::colnames(sVals$wald) = varNames
    base::rownames(sVals$wald) = base::paste0("Bootstrap",1:B)
    waldTime = tictoc::toc(quiet = TRUE)
    timing3[2] = waldTime$toc - waldTime$tic
  }
  if(dev == TRUE){
    tictoc::tic()
    devQD = base::matrix(data = NA,
                         nrow = B,
                         ncol = p)
    devPath = lapply(X = sValues,
                     FUN = solution_path,
                     row = 3)
    devSpace = reduced_space(list = devPath, vector = varNames)
    output$devSpace = devSpace
    sVals$dev = base::matrix(data = NA_integer_,
                             nrow = B,
                             ncol = p)
    for (i in 1:B){
      sVals$dev[i,] = sValues[[i]][[1]][3,]
      devQD[i,] = sValues[[i]][[2]]
    }
    base::colnames(sVals$dev) = varNames
    base::rownames(sVals$dev) = base::paste0("Bootstrap",1:B)
    devTime = tictoc::toc(quiet = TRUE)
    timing3[3] = devTime$toc - devTime$tic
  }

  countModels = count_models(output,
                             k = k)
  output = check_model_space(output,
                             k = k,
                             data = data,
                             family = family,
                             tcc = tcc)

  timing2 = base::rep(0,
                      times = 3)
  for(i in 1:base::length(output)){
    timing2[i] = output[[i]]$timing
    output[[i]] = output[[i]]$candidateSpace
  }
  base::names(output) = c("coef", "wald", "dev")[c(coef,wald,dev)]

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
  base::colnames(devQD) = varNames
  base::rownames(devQD) = base::paste0("B",
                                       1:B)
  output$devQD = devQD

  base::rownames(bootstraps) = base::paste0("B",
                                            1:B)
  output$bootsrap = bootstraps
  output$varNames = varNames

  finalTime = base::rep(0,
                        times = 3)
  for(b in 1:B){
    finalTime = finalTime + sValues[[b]][[3]]
  }

  finalTime = finalTime + timing2 + timing3 + timing

  output$time = finalTime

  return(output)
}

#####

#' Reduced model space
#'
#' @importFrom magrittr %>%



reduced_space = function(list,
                         vector){
  allModels = list %>%
    plyr::ldply(., data.frame) %>%
    plyr::ddply(., c("Variables", "Dimension"), dplyr::summarise,
                Count = length(Variables)) %>%
    dplyr::arrange(Dimension) %>%
    dplyr::mutate(Variables = paste(Variables, sep = "+"))
  return(allModels)
}

#####

#' @importFrom magrittr %>%

select_k = function(df,
                    k){
  df$Count = df$Count %>%
    as.numeric()
  df %>%
    dplyr::filter(Count >= k)
}

#####

#' @importFrom magrittr %>%

counting = function(list){
  test = matrix(data = 0, nrow  =3, ncol = 7)
  for(i in 1:3){
    check = list[[i]]%>% group_by(., dimension) %>% count()
    test[i,] = check$n
  }
  return(test)
}

#####

count_models = function(list,
                        k = 1){
  list = base::lapply(X = list,
                      FUN = select_k,
                      k = k)
  lapply(X = list,
         FUN = base::NROW)
}

#####

check_model_space = function(list,
                             k = 1,
                             data,
                             family,
                             tcc = tcc){
  list = base::lapply(X = list,
                      FUN = select_k,
                      k = k)
  output = base::lapply(X = list,
                        FUN = matrix_IC,
                        data = data,
                        family = family,
                        tcc = tcc)
  for( i in 1:length(output)){
    output[[i]]$candidateSpace = plyr::ldply (output[[i]]$candidateSpace,
                               base::data.frame,
                               stringsAsFactors = FALSE)
    rownames(output[[i]]$candidateSpace) = paste0("Model ",1:NROW(output[[i]]$candidateSpace))
  }
  #  output = base::lapply(X = output,
  #                        FUN = best_model)
  return(output)
}
