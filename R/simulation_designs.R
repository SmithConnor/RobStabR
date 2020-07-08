paper_results_simple = function(K = 100,
                                startSeed = 1234,
                                n,
                                p,
                                beta,
                                name,
                                family = poisson(link = "log")){
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
    tic()
    RobStab[[i]] = model_space(data = data,
                            B = 100,
                            m = n/2,
                            nStrata = 8,
                            family = family,
                            k = 1,
                            resid = "pearson",
                            coef = TRUE,
                            wald = TRUE,
                            dev = TRUE)
    t1 = toc()
    time[i,1] = t1$toc - t1$tic
    runs$RobStab = RobStab[[i]]
    tic()
    stepGLM[[i]] = step_glmrob(data = data,
                               family = family)
    t2 = toc()
    time[i,2] = t2$toc - t2$tic
    runs$stepQD = stepGLM[[i]]
    tic()
    stepAIC[[i]] = step_ic(data = data,
                           family = family,
                           k = 2)
    t3 = toc()
    time[i,3] = t3$toc - t3$tic
    runs$stepAIC = stepAIC[[i]]
    tic()
    stepBIC[[i]] = step_ic(data = data,
                           family = family,
                           k = log(n))
    t4 = toc()
    time[i,4] = t4$toc - t4$tic
    runs$stepBIC = stepBIC[[i]]
    tic()
    bestAIC[[i]] = bestglm_ic(data = data,
                              family = family,
                              IC = "AIC")
    t5 = toc()
    time[i,5] = t5$toc - t5$tic
    runs$bestAIC = bestAIC[[i]]
    tic()
    bestBIC[[i]] = bestglm_ic(data = data,
                              family = family,
                              IC = "BIC")
    t6 = toc()
    time[i,6] = t6$toc - t6$tic
    runs$bestBIC = bestBIC[[i]]
    tic()
    exhaustRDBC[[i]] = exhaustive_RDBC(data = data,
                                       family = family)
    t7 = toc()
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