gmp <- read.table("../data/gmp.dat")
gmp$pop <- round(gmp$gmp/gmp$pcgmp)

mse <- function(parameter_vector = vector(length = 2), N = gmp$pop, Y = gmp$pcgmp){
  y0 <- parameter_vector[1]
  a <- parameter_vector[2]
  mse <- mean((y0 * N^a - Y)^2)
  return(mse)
}

plm <- function(initial_y0, initial_a, N = gmp$pop, Y = gmp$pcgmp){
  result <- nlm(function(parameter_vector)mse(parameter_vector = vector(length = 2), N, Y), c(y0 = initial_y0, a = initial_a))
  final <- list("final_guess_y0" = result$estimate[1], "final_guess_a" = result$estimate[2], "final_value_MSE" = result$minimum)
  return(final)
}

plm.jackknife <- function(initial_y0, initial_a, N = gmp$pop, Y = gmp$pcgmp){
  n <- nrow(gmp)
  jackknife_parameter <- matrix(0, nrow = 2, ncol = n)
  for(omitted_point in 1:n){
    jackknife_parameter[1,omitted_point] <- plm(initial_y0, initial_a, N[-omitted_point], Y[-omitted_point])$final_guess_y0
    jackknife_parameter[2,omitted_point] <- plm(initial_y0, initial_a, N[-omitted_point], Y[-omitted_point])$final_guess_a
  }
  jackknife_var <- ((n-1)^2 / n) * apply(jackknife_parameter, 1, var)
  
  return(sqrt(jackknife_var))
}

suppressWarnings(plm.jackknife(6611, 1/8))