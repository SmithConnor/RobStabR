paper_results_simple = function(K = 100,
                                startSeed = 1234,
                                n,
                                p,
                                beta,
                                family = poisson(link = "log"),
                                tcc){
  yFinal = base::matrix(data = NA_integer_,
                        ncol = K,
                        nrow = n)
  seedSample = base::rep(x = NA_integer_,
                         times = K)
  poisSamp = base::list()
  set.seed(startSeed)
  x = mvtnorm::rmvnorm(n = n,
                       mean = base::rep(x = 0,
                                        times = p),
                       sigma = base::diag(p))
  eta = x %*% beta
  mu = exp(eta)
  RobStab = list()
  stepGLM = list()
  stepAIC = list()
  stepBIC = list()
  bestAIC = list()
  bestBIC = list()
  exhaustRDBC = list()
  time = matrix(data = 0, nrow = K, ncol = 7)
  for (i in 1:K){
    runs = list()
    seed = base::sample(x = 1:100000,
                        size = 1)
    print(seed)
    seedSample[i] = seed
    y = stats::rpois(n = n, lambda = mu)
    yFinal[,i] = y
    data = data.frame(x,y)
    tictoc::tic()
    RobStab[[i]] = model_space(data = data,
                            B = 100,
                            m = n/2,
                            nStrata = 8,
                            family = family,
                            k = 1,
                            resid = "pearson",
                            coef = TRUE,
                            wald = TRUE,
                            dev = TRUE,
                            tcc = tcc)
    t1 = tictoc::toc(quiet = TRUE)
    time[i,1] = t1$toc - t1$tic
    runs$RobStab = RobStab[[i]]
    tictoc::tic()
    stepGLM[[i]] = step_glmrob(data = data,
                               family = family,
                               tcc = tcc)
    t2 = tictoc::toc(quiet = TRUE)
    time[i,2] = t2$toc - t2$tic
    runs$stepQD = stepGLM[[i]]
    tictoc::tic()
    stepAIC[[i]] = step_ic(data = data,
                           family = family,
                           k = 2)
    t3 = tictoc::toc(quiet = TRUE)
    time[i,3] = t3$toc - t3$tic
    runs$stepAIC = stepAIC[[i]]
    tictoc::tic()
    stepBIC[[i]] = step_ic(data = data,
                           family = family,
                           k = log(n))
    t4 = tictoc::toc(quiet = TRUE)
    time[i,4] = t4$toc - t4$tic
    runs$stepBIC = stepBIC[[i]]
    tictoc::tic()
    bestAIC[[i]] = bestglm_ic(data = data,
                              family = family,
                              IC = "AIC")
    t5 = tictoc::toc(quiet = TRUE)
    time[i,5] = t5$toc - t5$tic
    runs$bestAIC = bestAIC[[i]]
    tictoc::tic()
    bestBIC[[i]] = bestglm_ic(data = data,
                              family = family,
                              IC = "BIC")
    t6 = tictoc::toc(quiet = TRUE)
    time[i,6] = t6$toc - t6$tic
    runs$bestBIC = bestBIC[[i]]
    tictoc::tic()
    exhaustRDBC[[i]] = exhaustive_RDBC(data = data,
                                       family = family,
                                       tcc = tcc)
    t7 = tictoc::toc(quiet = TRUE)
    time[i,7] = t7$toc - t7$tic
    runs$exhaustRDBC = exhaustRDBC[[i]]
    runs$time = time[i,]
  }
  results = base::list(RobStab = RobStab,
                 stepQD = stepGLM,
                 stepAIC = stepAIC,
                 stepBIC = stepBIC,
                 bestAIC = bestAIC,
                 bestBIC = bestBIC,
                 exhaustRDBC = exhaustRDBC,
                 x = x,
                 y  = y,
                 seed = seedSample,
                 beta = beta,
                 time = time)
  return(results)
}

############


paper_results_complex = function(K = 100, startSeed = 1234, n, p, beta, family = poisson(link = "log"), sigma = 0.2, sigma2 = 0.5, outlier, outlier2, k = 1, tcc){
  yFinal = base::matrix(data = NA_integer_,
                        ncol = K,
                        nrow = n)
  seedSample = base::rep(x = NA_integer_,
                         times = K)
  poisSamp = base::list()
  corMat = sigma^(toeplitz(1:p - 1))
  for(i in 1:p){
    for(j in 1:p){
      if( i %% 2 == 1 & j - i == 1 | j %% 2 == 1 & i - j == 1){
        corMat[i,j] = sigma2
      }
    }
  }
  set.seed(startSeed)
  x = mvtnorm::rmvnorm(n = n,
                       mean = base::rep(x = 0,
                                        times = p),
                       sigma = corMat)
  eta = x %*% beta
  mu = exp(eta)
  RSSM = list()
  stepGLM = list()
  stepAIC = list()
  stepBIC = list()
  time = matrix(data = 0, nrow = K, ncol = 4)
  for (i in 1:K){
    runs = list()
    seed = base::sample(x = 1:100000,
                        size = 1)
    seedSample[i] = seed
    mu[which(n - rank(x[,8]) < 60)] = outlier
    mu[which(n - rank(x[,16]) < 10)] = outlier2
    y = stats::rpois(n = n, lambda = mu)
    yFinal[,i] = y
    data = data.frame(x,y)
    tictoc::tic()
    RSSM[[i]] = model_space(data = data,
                            B = 100,
                            m = n/2,
                            nStrata = 8,
                            family = family,
                            k = k,
                            resid = "pearson",
                            coef = TRUE,
                            wald = TRUE,
                            dev = FALSE,
                            tcc = 2)
    t1 = tictoc::toc(quiet = TRUE)
    time[i,1] = t1$toc - t1$tic
    runs$RSSM = RSSM[[i]]
    tictoc::tic()
    stepGLM[[i]] = step_glmrob(data = data,
                               family = family,
                               tcc = 2)
    t2 = tictoc::toc(quiet = TRUE)
    time[i,2] = t2$toc - t2$tic
    runs$stepQD = stepGLM[[i]]
    tictoc::tic()
    stepAIC[[i]] = step_ic(data = data,
                           family = family,
                           k = 2)
    t3 = tictoc::toc(quiet = TRUE)
    time[i,3] = t3$toc - t3$tic
    runs$stepAIC = stepAIC[[i]]
    tictoc::tic()
    stepBIC[[i]] = step_ic(data = data,
                           family = family,
                           k = log(n))
    t4 = tictoc::toc(quiet = TRUE)
    time[i,4] = t4$toc - t4$tic
    runs$stepBIC = stepBIC[[i]]
    runs$time = time[i,]
  }
  results = list(RSSM = RSSM,
                 stepQD = stepGLM,
                 stepAIC = stepAIC,
                 stepBIC = stepBIC,
                 x = x,
                 y  = y,
                 seed = seedSample,
                 beta = beta,
                 time = time)
  return(results)
}



#####
paper_results_medium = function(K = 100, seedStart = 1234, n, p, beta, family = poisson(link = "log"), sigma = 0.2, outlier, tcc = 2){
  yFinal = base::matrix(data = NA_integer_,
                        ncol = K,
                        nrow = n)
  seedSample = base::rep(x = NA_integer_,
                         times = K)
  poisSamp = base::list()
  set.seed(seedStart)
  x = mvtnorm::rmvnorm(n = n,
                       mean = base::rep(x = 0,
                                        times = p),
                       sigma = sigma^(toeplitz(1:p - 1)))
  eta = x %*% beta
  mu = exp(eta)
  RSSM = list()
  stepGLM = list()
  stepAIC = list()
  stepBIC = list()
  bestAIC = list()
  bestBIC = list()
  exhaustRDBC = list()
  time = matrix(data = 0, nrow = K, ncol = 7)

  for (i in 1:K){
    runs = list()
    seed = base::sample(x = 1:100000,
                        size = 1)
    seedSample[i] = seed
    mu[which(n - rank(x[,5]) < 30)] = outlier
    y = stats::rpois(n = n, lambda = mu)
    yFinal[,i] = y
    data = data.frame(x,y)
    tictoc::tic()
    RSSM[[i]] = model_space(data = data,
                            B = 100,
                            m = n/2,
                            nStrata = 8,
                            family = family,
                            k = 1,
                            resid = "pearson",
                            coef = TRUE,
                            wald = TRUE,
                            dev = FALSE,
                            tcc = tcc)
    t1 = tictoc::toc(quiet = TRUE)
    time[i,1] = t1$toc - t1$tic
    runs$RSSM = RSSM[[i]]
    tictoc::tic()
    stepGLM[[i]] = step_glmrob(data = data,
                               family = family,
                               tcc = tcc)
    t2 = tictoc::toc(quiet = TRUE)
    time[i,2] = t2$toc - t2$tic
    runs$stepQD = stepGLM[[i]]
    tictoc::tic()
    stepAIC[[i]] = step_ic(data = data,
                           family = family,
                           k = 2)
    t3 = tictoc::toc(quiet = TRUE)
    time[i,3] = t3$toc - t3$tic
    runs$stepAIC = stepAIC[[i]]
    tictoc::tic()
    stepBIC[[i]] = step_ic(data = data,
                           family = family,
                           k = log(n))
    t4 = tictoc::toc(quiet = TRUE)
    time[i,4] = t4$toc - t4$tic
    runs$stepBIC = stepBIC[[i]]
    tictoc::tic()
    bestAIC[[i]] = bestglm_ic(data = data,
                              family = family,
                              IC = "AIC")
    t5 = tictoc::toc(quiet = TRUE)
    time[i,5] = t5$toc - t5$tic
    runs$bestAIC = bestAIC[[i]]
    tictoc::tic()
    bestBIC[[i]] = bestglm_ic(data = data,
                              family = family,
                              IC = "BIC")
    t6 = tictoc::toc(quiet = TRUE)
    time[i,6] = t6$toc - t6$tic
    runs$bestBIC = bestBIC[[i]]
    # tictoc::tic()
    # exhaustRDBC[[i]] = exhaustive_RDBC(data = data,
    #                                    family = family)
    # t7 = tictoc::toc(quiet = TRUE)
    # time[i,7] = t7$toc - t7$tic
    # runs$exhaustRDBC = exhaustRDBC[[i]]
    runs$time = time[i,]
  }
  results = list(RSSM = RSSM,
                 stepQD = stepGLM,
                 stepAIC = stepAIC,
                 stepBIC = stepBIC,
                 bestAIC = bestAIC,
                 bestBIC = bestBIC,
                 x = x,
                 y  = y,
                 seed = seedSample,
                 beta = beta,
                 time = time)
  return(results)
}


currentRobStab = example_3
beta = currentRobStab$beta
trueModel = paste0("X",which(beta != 0))
stepGLMRes = lapply(currentRobStab$stepQD, model_check, true = trueModel) %>%
  unlist() %>%
  matrix(., ncol = 3, byrow = TRUE)
stepAICRes = lapply(currentRobStab$stepAIC, model_check, true = trueModel) %>%
  unlist() %>%
  matrix(., ncol = 3, byrow = TRUE)
stepBICRes = lapply(currentRobStab$stepBIC, model_check, true = trueModel) %>%
  unlist() %>%
  matrix(., ncol = 3, byrow = TRUE)
bestAICRes = lapply(currentRobStab$bestAIC, model_check, true = trueModel) %>%
  unlist() %>%
  matrix(., ncol = 3, byrow = TRUE)
bestBICRes = lapply(currentRobStab$bestBIC, model_check, true = trueModel) %>%
  unlist() %>%
  matrix(., ncol = 3, byrow = TRUE)
exhaustRDBCRes = lapply(currentRobStab$exhaustRDBC, model_check, true = trueModel) %>%
  unlist() %>%
  matrix(., ncol = 3, byrow = TRUE)
RobStabP1 = lapply(currentRobStab$RobStab, best_robust, penalty = "P1") %>%
  lapply(., model_check_RobStab, true = trueModel)
matrixP1 = matrix(data = 0, ncol = 3, nrow = 3)
colnames(matrixP1) = c("wrong", "true", "correct")
rownames(matrixP1) = c("coef", "wald", "dev")
for( i in 1:K){
  matrixP1 = matrixP1 + RobStabP1[[i]]
}
RobStabP2 = lapply(currentRobStab$RobStab, best_robust, penalty = "P2") %>%
  lapply(., model_check_RobStab, true = trueModel)
matrixP2 = matrix(data = 0, ncol = 3, nrow = 3)
colnames(matrixP2) = c("wrong", "true", "correct")
rownames(matrixP2) = c("coef", "wald", "dev")
for( i in 1:K){
  matrixP2 = matrixP2 + RobStabP2[[i]]
}

matrixP1
matrixP2
stepGLMRes %>% apply(.,2,sum)
stepAICRes %>% apply(.,2,sum)
stepBICRes %>% apply(.,2,sum)
bestAICRes %>% apply(.,2,sum)
bestBICRes %>% apply(.,2,sum)
exhaustRDBCRes %>% apply(.,2,sum)
#
