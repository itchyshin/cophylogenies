# coustm functions

#' Title: getting Zr and its sampling variance from p value 
#'
#' @param data: data frame 
#' @param pval: p value
#' @param N: sample size (N: the number of species ) and the degrees of freedom df = N - 2
#'
#' @return
#' @export
#'
#' @examples
p_to_Zr <- function(data, pval, N) {
  
  # turning them into strings
  pval <- data[[deparse(substitute(pval))]]
  N <- data[[deparse(substitute(N))]]
  
  # getting t values 
  tval<- -qt(pval, N - 2) 
  rval <- tval / sqrt((tval^2) + (N - 2))
  
  # define Zr function
  # Zr <- 0.5*(log(1 + rval) - log(1 - rval)); the same as below
  # r <-tanh(Zr) # turning Zr to r
  Zr <- atanh(rval)
  
  # getting Var(Zr)
  VZr <- 1 / (N - 3)
  
  # putting all together
  Zrs <- tibble(rval, Zr, VZr)
  data <- bind_cols(data, Zrs)
}

# coverting back Zr to r
# Just use "psych" pacakge - fisherz2r(z) - <http://personality-project.org/r/psych/help/fisherz.html>
# or this will do : r to Zr is tanh(r)!!

# Functions for processing


# General modeling functions 
# Functions for I2

#' Title Function to obtain total and separate I2 from multilevel-meta-analytic model
#'
#' @param model 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
I2 <- function(model, method = c("Wolfgang", "Shinichi")){
  
  ## evaluate choices
  method <- match.arg(method)
  
  # Wolfgang's method
  if(method == "Wolfgang"){
    W <- solve(model$V) 
    X <- model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    I2_each  <- model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
    
    # or my way
  } else {
    # sigma2_v = typical sampling error variance
    sigma2_v <- sum(1/model$vi) * (model$k-1) / (sum(1/model$vi)^2 - sum((1/model$vi)^2)) 
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + sigma2_v) #s^2_t = total variance
    I2_each  <- model$sigma2 / (sum(model$sigma2) + sigma2_v)
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
  }
  return(I2s)
}

# test <- dataset$fit4.1[[3]]
# I2(test, method = "Wolfgang")
# I2(test, method = "Shinichi")


#' Title: R2 based on Nakagawa & Schielzeth 2013
#'
#' @param model 
#'
#' @return
#' @export
#'
#' @examples
R2 <- function(model){
  warning("Conditional R2 is not meaningful and the same as marginal R2\n")
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  
  # marginal
  R2m <- fix / (fix + sum(model$sigma2))
  R2
  #Rm <- round(100*R2m, 3)
  
  # conditional
  R2c <- (fix + sum(model$sigma2) - model$sigma2[length(model$sigma2)]) / 
    (fix + sum(model$sigma2))
  
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  return(R2s)
}


#' Title: the function to get estimates from rma objects (metafor)
#'
#' @param model: rma.mv object 
#' @param mod: the name of a moderator 
get_est <- function (model, mod = " ") {
  
  name <- as.factor(str_replace(row.names(model$beta), mod, ""))
  estimate <- as.numeric(model$beta)
  lowerCL <- model$ci.lb
  upperCL <- model$ci.ub 
  
  table <- tibble(name = name, estimate = estimate, lowerCL = lowerCL, upperCL = upperCL)
}


#' Title: the function to get prediction intervals (crediblity intervals) from rma objects (metafor)
#'
#' @param model: rma.mv object 
#' @param mod: the name of a moderator 
get_pred <- function (model, mod = " ") {
  name <- as.factor(str_replace(row.names(model$beta), mod, ""))
  len <- length(name)
  
  if(len != 1){
  newdata <- matrix(NA, ncol = len, nrow = len)
  for(i in 1:len) {
    # getting the position of unique case from X (design matrix)
    pos <- which(model$X[,i] == 1)[[1]]
    newdata[, i] <- model$X[pos,]
    }
  pred <- predict.rma(model, newmods = newdata)
  }
  else {
    pred <- predict.rma(model)
    }
  lowerPR <- pred$cr.lb
  upperPR <- pred$cr.ub 
  
  table <- tibble(name = name, lowerPR = lowerPR, upperPR = upperPR)
}

#Here are links for how to do confidence regions for rma.mv regression lines
#https://www.rdocumentation.org/packages/metafor/versions/1.9-9/topics/predict.rma
#https://stackoverflow.com/questions/50804464/out-of-sample-prediction-for-rma-object-in-metafor


#' Title: Contrast name geneator
#'
#' @param name: a vector of character strings
cont_gen <- function (name) {
  combination <- combn(name,2)
  name_dat <- t(combination)
  names <- paste(name_dat[ ,1], name_dat[, 2], sep = "-")
  return(names)
}
