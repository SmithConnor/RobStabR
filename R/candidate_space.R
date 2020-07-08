#' Find the candidate model space after B bootstraps
#'
#' @param data a data frame containing the variables in the model.
#' @param B
#' @param m
#' @param nStrata
#' @param family
#' @param k
#' @param resid
#' @param coef
#' @param wald
#' @param dev
#' @param bootstraps
#' @importFrom magrittr %>%
#'
#' @example



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
                       bootstraps = NA){
  # Parameters
  n = base::NROW(data)
  p = base::NCOL(data) - 1
  nMethods = coef + wald + dev
  # Variables
  devQD = 0
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
                               control = robustbase::glmrobMqle.control(tcc = 2,
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
                        .parallel = FALSE)
  if(coef == TRUE){
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
  }
  if(wald == TRUE){
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
  }
  if(dev == TRUE){
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
  }

  countModels = count_models(output,
                             k = k)
  output = check_model_space(output,
                             k = k,
                             data = data,
                             family = family)
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
  return(output)
}

#####

#' Reduced model space
#'
#' @param list
#' @param vector
#' @importFrom magrittr %>%
#'
#' @example


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

#' @param df
#' @param k
#' @importFrom magrittr %>%

select_k = function(df,
                    k){
  df$Count = df$Count %>%
    as.numeric()
  df %>%
    dplyr::filter(Count >= k)
}

#####

#' @param list
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

#' @param data a data frame containing the variables in the model.

check_model_space = function(list,
                             k = 1,
                             data,
                             family){
  list = base::lapply(X = list,
                      FUN = select_k,
                      k = k)
  output = base::lapply(X = list,
                        FUN = matrix_IC,
                        data = data,
                        family = family)
  for( i in 1:length(output)){
    output[[i]] = plyr::ldply (output[[i]],
                               base::data.frame,
                               stringsAsFactors = FALSE)
    rownames(output[[i]]) = paste0("Model ",1:NROW(output[[i]]))
  }
  #  output = base::lapply(X = output,
  #                        FUN = best_model)
  return(output)
}
