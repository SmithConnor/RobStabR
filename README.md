
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
    ## Model 1                    1+X1         1     4 361.36265 362.36265
    ## Model 2                    1+X2         1     2 289.41846 290.41846
    ## Model 3                    1+X3         1    94 115.36365 116.36365
    ## Model 4                 1+X1+X3         2    38 116.56831 118.56831
    ## Model 5                 1+X2+X3         2    62  40.44998  42.44998
    ## Model 6              1+X1+X2+X3         3    74  19.31966  22.31966
    ## Model 7              1+X1+X3+X4         3     4 114.64998 117.64998
    ## Model 8              1+X1+X3+X5         3     3 120.32996 123.32996
    ## Model 9              1+X1+X3+X6         3     1 120.66085 123.66085
    ## Model 10             1+X1+X3+X7         3     3 111.48272 114.48272
    ## Model 11             1+X1+X3+X8         3     2 111.49674 114.49674
    ## Model 12             1+X2+X3+X4         3     1  44.53229  47.53229
    ## Model 13             1+X2+X3+X5         3     1  44.34436  47.34436
    ## Model 14             1+X2+X3+X6         3     3  44.05728  47.05728
    ## Model 15             1+X2+X3+X7         3     8  43.35645  46.35645
    ## Model 16          1+X1+X2+X3+X4         4    25  22.06696  26.06696
    ## Model 17          1+X1+X2+X3+X5         4    14  23.32576  27.32576
    ## Model 18          1+X1+X2+X3+X6         4    15  22.72866  26.72866
    ## Model 19          1+X1+X2+X3+X7         4    19  20.61659  24.61659
    ## Model 20          1+X1+X2+X3+X8         4    17  23.42801  27.42801
    ## Model 21          1+X1+X3+X4+X7         4     2 108.62182 112.62182
    ## Model 22          1+X1+X3+X5+X7         4     2 113.87718 117.87718
    ## Model 23          1+X1+X3+X5+X8         4     1 113.21072 117.21072
    ## Model 24          1+X2+X3+X6+X7         4     2  46.82274  50.82274
    ## Model 25          1+X2+X3+X7+X8         4     3  46.77629  50.77629
    ## Model 26       1+X1+X2+X3+X4+X5         5    19  25.55239  30.55239
    ## Model 27       1+X1+X2+X3+X4+X6         5     6  25.73329  30.73329
    ## Model 28       1+X1+X2+X3+X4+X7         5     7  23.02295  28.02295
    ## Model 29       1+X1+X2+X3+X4+X8         5    13  26.07085  31.07085
    ## Model 30       1+X1+X2+X3+X5+X6         5     6  26.78906  31.78906
    ## Model 31       1+X1+X2+X3+X5+X7         5    12  24.15976  29.15976
    ## Model 32       1+X1+X2+X3+X5+X8         5     9  27.35279  32.35279
    ## Model 33       1+X1+X2+X3+X6+X7         5     8  23.75736  28.75736
    ## Model 34       1+X1+X2+X3+X6+X8         5     8  26.84586  31.84586
    ## Model 35       1+X1+X2+X3+X7+X8         5     9  24.61842  29.61842
    ## Model 36       1+X1+X3+X4+X5+X7         5     2 107.15546 112.15546
    ## Model 37       1+X2+X3+X4+X6+X7         5     1  50.86694  55.86694
    ## Model 38    1+X1+X2+X3+X4+X5+X6         6     8  29.36426  35.36426
    ## Model 39    1+X1+X2+X3+X4+X5+X7         6     8  25.40305  31.40305
    ## Model 40    1+X1+X2+X3+X4+X5+X8         6    21  29.14314  35.14314
    ## Model 41    1+X1+X2+X3+X4+X6+X7         6    18  26.48768  32.48768
    ## Model 42    1+X1+X2+X3+X4+X6+X8         6     3  29.76000  35.76000
    ## Model 43    1+X1+X2+X3+X4+X7+X8         6    11  27.12182  33.12182
    ## Model 44    1+X1+X2+X3+X5+X6+X7         6    12  27.40829  33.40829
    ## Model 45    1+X1+X2+X3+X5+X6+X8         6     7  30.84692  36.84692
    ## Model 46    1+X1+X2+X3+X5+X7+X8         6     4  28.27613  34.27613
    ## Model 47    1+X1+X2+X3+X6+X7+X8         6     7  27.71395  33.71395
    ## Model 48    1+X1+X3+X4+X5+X7+X8         6     1 102.39641 108.39641
    ## Model 49 1+X1+X2+X3+X4+X5+X6+X7         7    23  29.12117  36.12117
    ## Model 50 1+X1+X2+X3+X4+X5+X6+X8         7    19  33.03842  40.03842
    ## Model 51 1+X1+X2+X3+X4+X5+X7+X8         7    31  29.53256  36.53256
    ## Model 52 1+X1+X2+X3+X4+X6+X7+X8         7    12  30.55453  37.55453
    ## Model 53 1+X1+X2+X3+X5+X6+X7+X8         7    15  31.49186  38.49186

The best model selected byt the RDBC when the penalty is P2 is:

``` r
robStabResult$dev[which.min(robStabResult$dev$P2),]
```

    ##              model dimension count       P1       P2
    ## Model 6 1+X1+X2+X3         3    74 19.31966 22.31966
