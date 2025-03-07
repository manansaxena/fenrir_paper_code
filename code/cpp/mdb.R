# this code is specifically for checking how multinomial dirichlet bootstrap compares to our Debiased Multinomial Dirichlet Bootstrap and Stan

library(MCMCpack)
library(compositions)

multinomial_dirichlet_bootstrap <- function(Y, alpha = 1) {

  D <- nrow(Y)
  N <- ncol(Y)
  at <- colSums(Y) + alpha * D
  E_pi_t <- sweep(Y + alpha, 2, at, "/")
  E_pi_t <- apply(E_pi_t, 2, function(p) alr(p))
  return(E_pi_t)
}

set.seed(123)
D <- 10
Q <- 1
data_path <- "random_path_where_data_is_stored"
result_path <- "random_path_where_results_are_stored"
n_cores <- 16

# Load data
data <- readRDS(paste0(data_path, "data.rds"))
Y_obs <- data$Y_obs_combined
N_obs <- dim(Y_obs)[2]
E_pi_t <- multinomial_dirichlet_bootstrap(Y_obs, alpha = 0.5)
write.csv(E_pi_t, file = paste0(result_path,"eta.csv"), row.names = FALSE)
