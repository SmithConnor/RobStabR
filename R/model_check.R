#' A function to check whether the selected model is a correct or wrong model.
#' @param model a vector containing which variables were selected.
#' @param true a vector containing all the variables in the true model.
#' @export

model_check = function(model,
                       true){
  vector = rep(NA, times = 3)
  dim = length(true)
  if(sum(true %in% model) == dim){
    if(length(model) == dim){
      vector = c(FALSE, TRUE, FALSE)
    } else{
      vector = c(FALSE, FALSE, TRUE)
    }
  } else{
    vector = c(TRUE, FALSE, FALSE)
  }
  names(vector) = c("wrong", "true", "correct")
  return(vector)
}

#####

#' A function to check whether the selected model is a correct or wrong model for the RobStab method.
#' @param list a list containing vectors for each subtractive lack-of-fit choice and which variables were selected.
#' @param true a vector containing all the variables in the true model.
#' @export

model_check_RobStab = function(list,
                               true){
  matrix = matrix(data = NA, ncol = 3, nrow = 3)
  colnames(matrix) = c("wrong", "true", "correct")
  rownames(matrix) = c("coef", "wald", "dev")
  for(j in 1:3){
    model = list[[j]]
    vector = rep(NA, times = 3)
    dim = length(true)
    if(sum(true %in% model) == dim){
      if(length(model) == dim){
        vector = c(FALSE, TRUE, FALSE)
      } else{
        vector = c(FALSE, FALSE, TRUE)
      }
    } else{
      vector = c(TRUE, FALSE, FALSE)
    }
    matrix[j,] = vector
  }
  return(matrix)
}

#####

check_size = function(list){
  vector = rep(NA, 3)
  names(vector) = c("coef", "wald", "dev")
  for(j in 1:3){
    vector[j] = nrow(list[[j]])
  }
  return(vector)
}
