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

#' Title
#'
#' @param model 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
I2 <- function(model, method = c("Wolfgang", "Shinichi")){
  warning("Make sure you have the observation (effec size) level random effect\n")
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

# Functions for R2
# make a version which provides both R types of R^2

#' Title
#'
#' @param model 
#'
#' @return
#' @export
#'
#' @examples
R2 <- function(model){
  warning("Make sure you have the observation (effec size) level random effect as the last in the formula\n")
  
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

# R2(test)

# functions for taking out estimates (b) and CIs along with the name - the 3 models at the same time (use facet_wrap?)

# function to get estimates etc from rma.mv objects 
# currently only works when you have one moderator - need to change it more general

#TODO - we should get this to estimate N at this stage for all the levels


#' Title: the function to get estimates
#'
#' @param model: rma.mv object 
#' @param mod: the name of a moderator 
get_est <- function (model, mod = " ") {
  
  name <- as.factor(str_replace(row.names(model$beta), mod, ""))
  estimate <- as.numeric(model$beta)
  lowerCL <- model$ci.lb
  upperCL <- model$ci.ub 
  
  tibble(name = name, esZr = estimate, lciZr  = lowerCL, uclZr = upperCL) %>% 
    mutate(esr = tanh(esZr), lcir = tanh(lciZr), ucir = tanh(uclZr)) -> table
}


# making forest plots

#' Title: the function to plot categorical data 
#'
#' @param data 
#' @param title 
#'
#' @return
#' @export
#'
#' @examples
sum_plot <- function(data, title = "") {
  ggplot(data, aes(x = name, y = estimate, colour = name)) +
    geom_pointrange(aes(ymin = lowerCL, ymax = upperCL), size = 1.5) +
    geom_hline(yintercept = 0, lty = 2, lwd = 1, alpha = 0.5) +
    labs(x = "", y = "Zr (effect size)", title = title) +
    #scale_y_continuous(limits = c(-0.2, 1)) + # need to be careful - points will not show if they are not in this range
    coord_flip() +
    facet_grid(mdata ~.) +
    theme(#panel.background = element_rect(fill = "white"),
      axis.line.x = element_line(size = 2, colour = "black"),
      axis.ticks.x = element_line(size = 2),
      axis.ticks.length = unit(0.3,"cm"),
      axis.text.x = element_text(size = 15),
      axis.ticks.y = element_line(size = 0, colour = "white"),
      axis.text.y = element_text(size = 15, face = "plain", colour = "black"),
      title = element_text(size = 15, face = "italic"),
      axis.title.x = element_text(size = 15, face = "plain"),
      legend.position = "none")
}

