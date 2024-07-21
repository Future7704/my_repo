# 读取数据
gmp <- read.table("data/gmp.dat", header = TRUE)
gmp$pop <- round(gmp$gmp/gmp$pcgmp)

# 定义MSE函数
mse <- function(parameter_vector, N = gmp$pop, Y = gmp$pcgmp){
  y0 <- parameter_vector[1]
  a <- parameter_vector[2]
  mse_value <- mean((y0 * N^a - Y)^2)
  return(mse_value)
}

# 定义plm函数
plm <- function(initial_y0, initial_a, N = gmp$pop, Y = gmp$pcgmp){
  result <- nlm(function(parameter_vector) mse(parameter_vector, N, Y), c(y0 = initial_y0, a = initial_a))
  final <- list("final_guess_y0" = result$estimate[1], "final_guess_a" = result$estimate[2], "final_value_MSE" = result$minimum)
  return(final)
}

# 定义plm.jackknife函数
plm.jackknife <- function(initial_y0, initial_a, N = gmp$pop, Y = gmp$pcgmp){
  n <- length(N)  # 使用N的长度来确定样本数量
  jackknife_parameter <- matrix(0, nrow = 2, ncol = n)
  for (omitted_point in 1:n){
    jackknife_parameter[1, omitted_point] <- plm(initial_y0, initial_a, N[-omitted_point], Y[-omitted_point])$final_guess_y0
    jackknife_parameter[2, omitted_point] <- plm(initial_y0, initial_a, N[-omitted_point], Y[-omitted_point])$final_guess_a
  }
  jackknife_var <- ((n - 1)^2 / n) * apply(jackknife_parameter, 1, var)
  return(sqrt(jackknife_var))
}

# 运行jackknife并抑制警告信息
suppressWarnings(plm.jackknife(6611, 1/8))

