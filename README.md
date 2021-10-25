
# RobStabR

RobStabR is a robust model selection framework for Generalised Linear
Models by Connor Smith, Boris Guennewig and Samuel Müller.

This GitHub is a repository for the R package RobStab R.

# Installation

``` r
library(devtools)
devtools::install_github("SmithConnor/RobStabR")
```

# An example

We have data with 64 observations and 8 variables, and we desire to
perform variables selection using the poisson model. We induce outliers
in the 8 observations where X8 is the largest by setting lambda = 5.

``` r
library(RobStabR)

set.seed(15)
n = 64
p = 8
beta = c(1, -1, 1, 
         rep(0, p - 3))
x = matrix(rnorm(n*p), 
           nrow = n,
           ncol = p)
colnames(x) = paste0("X", 
                     1:p)
lambda = exp(x %*% beta)
lambda[which(rank(x[,8]) <= 8)] = 5 # Induce outliers
y = rpois(n = n, 
          lambda = lambda)
data = data.frame(x, y)
robStabResult = RobStabR::model_space(data = data,
                                      B = 100,
                                      m = n/2,
                                      family = poisson(link = "log"),
                                      tcc = 2)

robStabResult$dev
```

    ##                           model dimension count        P1        P2
    ## Model 1                    1+X1         1     4 361.36892 362.36892
    ## Model 2                    1+X2         1     2 289.42180 290.42180
    ## Model 3                    1+X3         1    94 115.35614 116.35614
    ## Model 4                 1+X1+X3         2    38 116.56091 118.56091
    ## Model 5                 1+X2+X3         2    62  40.44450  42.44450
    ## Model 6              1+X1+X2+X3         3    74  19.31762  22.31762
    ## Model 7              1+X1+X3+X4         3     4 114.64133 117.64133
    ## Model 8              1+X1+X3+X5         3     3 120.32273 123.32273
    ## Model 9              1+X1+X3+X6         3     1 120.65342 123.65342
    ## Model 10             1+X1+X3+X7         3     3 111.47798 114.47798
    ## Model 11             1+X1+X3+X8         3     2 111.49064 114.49064
    ## Model 12             1+X2+X3+X4         3     1  44.52701  47.52701
    ## Model 13             1+X2+X3+X5         3     1  44.33877  47.33877
    ## Model 14             1+X2+X3+X6         3     3  44.05190  47.05190
    ## Model 15             1+X2+X3+X7         3     8  43.35179  46.35179
    ## Model 16          1+X1+X2+X3+X4         4    25  22.06403  26.06403
    ## Model 17          1+X1+X2+X3+X5         4    14  23.32395  27.32395
    ## Model 18          1+X1+X2+X3+X6         4    15  22.72688  26.72688
    ## Model 19          1+X1+X2+X3+X7         4    19  20.61629  24.61629
    ## Model 20          1+X1+X2+X3+X8         4    17  23.42609  27.42609
    ## Model 21          1+X1+X3+X4+X7         4     2 108.61558 112.61558
    ## Model 22          1+X1+X3+X5+X7         4     2 113.87333 117.87333
    ## Model 23          1+X1+X3+X5+X8         4     1 113.20559 117.20559
    ## Model 24          1+X2+X3+X6+X7         4     2  46.81821  50.81821
    ## Model 25          1+X2+X3+X7+X8         4     3  46.77148  50.77148
    ## Model 26       1+X1+X2+X3+X4+X5         5    19  25.54974  30.54974
    ## Model 27       1+X1+X2+X3+X4+X6         5     6  25.73065  30.73065
    ## Model 28       1+X1+X2+X3+X4+X7         5     7  23.02167  28.02167
    ## Model 29       1+X1+X2+X3+X4+X8         5    13  26.06808  31.06808
    ## Model 30       1+X1+X2+X3+X5+X6         5     6  26.78748  31.78748
    ## Model 31       1+X1+X2+X3+X5+X7         5    12  24.16019  29.16019
    ## Model 32       1+X1+X2+X3+X5+X8         5     9  27.35125  32.35125
    ## Model 33       1+X1+X2+X3+X6+X7         5     8  23.75743  28.75743
    ## Model 34       1+X1+X2+X3+X6+X8         5     8  26.84420  31.84420
    ## Model 35       1+X1+X2+X3+X7+X8         5     9  24.61805  29.61805
    ## Model 36       1+X1+X3+X4+X5+X7         5     2 107.15044 112.15044
    ## Model 37       1+X2+X3+X4+X6+X7         5     1  50.86269  55.86269
    ## Model 38    1+X1+X2+X3+X4+X5+X6         6     8  29.36185  35.36185
    ## Model 39    1+X1+X2+X3+X4+X5+X7         6     8  25.40271  31.40271
    ## Model 40    1+X1+X2+X3+X4+X5+X8         6    21  29.14097  35.14097
    ## Model 41    1+X1+X2+X3+X4+X6+X7         6    18  26.48680  32.48680
    ## Model 42    1+X1+X2+X3+X4+X6+X8         6     3  29.75751  35.75751
    ## Model 43    1+X1+X2+X3+X4+X7+X8         6    11  27.12051  33.12051
    ## Model 44    1+X1+X2+X3+X5+X6+X7         6    12  27.40905  33.40905
    ## Model 45    1+X1+X2+X3+X5+X6+X8         6     7  30.84559  36.84559
    ## Model 46    1+X1+X2+X3+X5+X7+X8         6     4  28.27646  34.27646
    ## Model 47    1+X1+X2+X3+X6+X7+X8         6     7  27.71393  33.71393
    ## Model 48    1+X1+X3+X4+X5+X7+X8         6     1 102.39271 108.39271
    ## Model 49 1+X1+X2+X3+X4+X5+X6+X7         7    23  29.12114  36.12114
    ## Model 50 1+X1+X2+X3+X4+X5+X6+X8         7    19  33.03644  40.03644
    ## Model 51 1+X1+X2+X3+X4+X5+X7+X8         7    31  29.53230  36.53230
    ## Model 52 1+X1+X2+X3+X4+X6+X7+X8         7    12  30.55361  37.55361
    ## Model 53 1+X1+X2+X3+X5+X6+X7+X8         7    15  31.49248  38.49248

The best model selected byt the RDBC when the penalty is P2 is:

``` r
robStabResult$dev[which.min(robStabResult$dev$P2),]
```

    ##              model dimension count       P1       P2
    ## Model 6 1+X1+X2+X3         3    74 19.31762 22.31762
