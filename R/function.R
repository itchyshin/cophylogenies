# function to calcuate a p value given correlation and DF


r_to_p <- function(r = 0.5, N = 23) {
  t_value <- sqrt( (r^2*(N - 2) ) / (1-r^2))
  p_value <- 1- abs(pt(t_value, N - 2))
  print(p_value)
}