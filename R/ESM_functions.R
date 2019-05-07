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

