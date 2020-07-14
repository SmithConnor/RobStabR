#' sMPP Plot
#'
#' @param robStab RobStabR object.
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot
#' @export

sMPP = function(robStab){
  output = list()
  for(i in 1:3){
    sVal = robStab[[i + 6]] %>%
      abs() %>%
      apply(.,1,rank)
    B = sVal %>%
      ncol()
    p = sVal %>%
      nrow()
    hMat = base::matrix(data = 0,
                        nrow = p,
                        ncol = p)
    for(j in 1:B){
      hMat = hMat + h_mat(sVal[,j])
    }
    colnames(hMat) = robStab$varNames
    rownames(hMat) = robStab$varNames
    ordering = sVal %>%
      apply(.,1,mean) %>%
      abs() %>%
      order(., decreasing = TRUE)
    hMat = hMat[ordering, ordering]
    hMat = hMat/B
    output[[i]] = hMat
  }
  base::names(output) = c("coef",
                          "wald",
                          "dev")
  output$plots = output %>%
    base::lapply(.,sMPPplot)

  return(output)
}

#####

h_mat = function(vector){
  p = vector %>%
    length()
  mat = base::matrix(data = 0,
                     nrow = p,
                     ncol = p)
  for(i in 1:p){
    for(j in 1:p){
      mat[i,j] = vector[i] > vector[j]
    }
  }
  return(mat)
}

###


sMPPplot = function(hMat){
  df = hMat %>%
    reshape2::melt()
  base::colnames(df) = c("Y", "X", "Probability")
  MPP = ggplot2::ggplot(data = df,
               mapping = ggplot2::aes(x = X,
                             y = reorder(Y, dplyr::desc(Y)))) +
    ggplot2::geom_tile(ggplot2::aes(fill = Probability), colour = "grey50") +
    ggplot2::scale_fill_gradient(low = "white", high = "black") +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::labs(y = "Variable Y", x = "Variable X") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face="bold",
                                                       color="#993333",
                                                       size = 16),
          axis.text.y = ggplot2::element_text(face="bold",
                                              color="#993333",
                                              size = 16),
          axis.text = ggplot2::element_text(face = "bold",
                                            size=16))
  return(MPP)
}

###

#' sVIP Plot
#' @param robStab RobStabR object.
#' @export
#'

sVIP = function(robStab){
  newRobStab = list(robStab[[7]],
                    robStab[[8]],
                    robStab[[9]])
  gMat = newRobStab %>%
    base::lapply(., g_mat)

  output = gMat

  base::names(output) = c("coef",
                          "wald",
                          "dev")

  output$plots = output %>%
    base::lapply(.,sVIPplot)

  return(output)
}

#####

sVIPplot = function(gMat){
  df = gMat %>%
    reshape2::melt()
  base::colnames(df) = c("Variable", "Dimension", "Probability")
  df$Dimension =  factor(df$Dimension, levels = base::colnames(gMat))
  base::levels(df$Dimension) = 0:base::nrow(gMat)
  VIP = ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = Dimension,
                                                          y = Probability,
                                                          col = Variable,
                                                          group = Variable)) +
    ggplot2::geom_line(lwd = 1.3)  +
    ggplot2::geom_point(lwd = 3) +
    ggplot2::labs(y = "Variable Inclusion Probability", x = "Model Dimension") +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face="bold",
                                                       color="#993333",
                                                       size = 16),
          axis.text.y = ggplot2::element_text(face="bold",
                                     color="#993333",
                                     size = 16),
          axis.text = ggplot2::element_text(face = "bold",
                                   size=16)) +
    ggplot2::scale_x_discrete(limits = base::factor(base::nrow(gMat):0))
  return(VIP)
}

#####

g_mat = function(matrix){
  p = base::NCOL(matrix)
  B = base::NROW(matrix)

  gMatrix = base::matrix(data = NA_real_,
                         nrow = p,
                         ncol = p + 1)

  colnames(gMatrix) = 0:p
  rownames(gMatrix) = colnames(matrix)

  sAbs = base::t(base::apply(X = matrix,
                             MARGIN = 1,
                             FUN = abs))
  sRank = p + 1 - base::t(apply(X = sAbs, MARGIN = 1, FUN = rank))

  for(i in 1:(p + 1)){
    for(j in 1:p){
      gMatrix[j, i] = base::sum(sRank[,j] <= (i - 1))
    }
  }

  gMatrix = gMatrix/B

  base::colnames(gMatrix) = base::paste0("Dim",0:p)
  return(gMatrix)
}
