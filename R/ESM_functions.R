# coustm functions

#' Title: getting Zr and its sampling variance from p value 
#'
#' @param data: data frame 
#' @param pval: p value
#' @param N: the number of species (N) and the degrees of freedom df = N - 2
#' @param n: alternative sample size: the number of species (n) and the degrees of freedom df = N - 2
#'
#' @return
#' @export
#'
#' @examples
p_to_Zr <- function(data, pval, N, n) {
  
  # turning them into strings
  pval <- data[[deparse(substitute(pval))]]
  N <- data[[deparse(substitute(N))]]
  n <- data[[deparse(substitute(n))]]
  
  # getting t values 
  tval<- -qt(pval, N - 2) 
  rval <- tval / sqrt((tval^2) + (N - 2))
  
  # define Zr function
  # Zr <- 0.5*(log(1 + rval) - log(1 - rval)); the same as below
  # r <-tanh(Zr) # turning Zr to r
  Zr <- atanh(rval)
  
  # getting Var(Zr)
  VZr <- 1 / (N - 3)
  
  # getting ZrA
  tvalA<- -qt(pval, n - 2) 
  rvalA <- tval / sqrt((tvalA^2) + (n - 2))
  ZrA <- 0.5*(log(1 + rvalA) - log(1 - rvalA))
  
  # Var(ZrA)
  VZrA <- 1 / (n - 3)
  
  # putting all together
  Zrs <- tibble(rval, Zr, VZr, rvalA, ZrA, VZrA)
  data <- bind_cols(data, Zrs)
}

# coverting back Zr to r
# Just use "psych" pacakge - fisherz2r(z) - <http://personality-project.org/r/psych/help/fisherz.html>
# or this will do : r to Zr is tanh(r)!!

