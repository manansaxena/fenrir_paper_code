# set the path where all the R libraries are stored.
.libPaths(c("random_path", .libPaths()))

library(fenrir)
library(fido)
library(coda)
library(abind)
library(driver)
library(compositions)
library(foreach)
library(doParallel)

optimizer <- function(Y_obs, observed_TT, N_total_list, F, G, gamma, W, M0, C0, Xi0, v0, init, result_path) {
  start_time <- Sys.time()
  res <- fenrir::fenrir_optim(
    Y_obs = Y_obs, observed_TT = observed_TT, N_total_list = N_total_list, F = F, G = G, gamma = gamma, W = W, M0 = M0, C0 = C0, Xi0 = Xi0, v0 = v0, init = init,
    log_probs_path = paste0(result_path, "logprobs_cpp.csv"), num_dirsamples = 2000, pseudocount = 0.5, max_iter = 10000, eps_f = 1e-5, eps_g = 1e-5
  )
  end_time <- Sys.time()
  res$total_time <- end_time - start_time
  return(res)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Please provide D, Q, data_path, result_path, and cores as command line arguments.")
}

# Parse the arguments

D <- as.integer(args[1])
Q <- as.integer(args[2])
data_path <- args[3]
result_path <- args[4]
n_cores <- as.integer(args[5])

# Load data
data <- readRDS(paste0(data_path, "data.rds"))

Y_obs <- data$Y_obs_combined
observed_TT <- data$observed_indices_combined

N_total <- length(observed_TT)
N_obs <- dim(Y_obs)[2]
P <- D - 1
num_timeseries <- length(data$N_total_list)
N_total_list <- data$N_total_list
M0_list <- data$M0_list
C0_list <- data$C0_list
W_val <- data$W_val
W <- lapply(1:N_total, function(x) diag(Q) * W_val)
F <- matrix(1, nrow = Q, ncol = N_total)
G <- lapply(1:N_total, function(x) diag(Q) * 1)
gamma <- rep(1, N_total)
M0 <- lapply(1:num_timeseries, function(x) matrix(0, nrow = Q, ncol = P))
replace_upper_tri <- function(mat, val) {
  mat[upper.tri(mat)] <- val
  return(mat)
}
M0 <- mapply(replace_upper_tri, M0, M0_list, SIMPLIFY = FALSE)
C0 <- lapply(1:num_timeseries, function(x) diag(Q))
C0 <- mapply(function(mat, val) {
  mat[] <- val
  mat
}, C0, C0_list, SIMPLIFY = FALSE)
Xi0 <- 1 * diag(P)
v0 <- D + 3
init <- matrix(0, nrow = P, ncol = N_obs)

result_optim <- optimizer(Y_obs, observed_TT, N_total_list, F, G, gamma, W, M0, C0, Xi0, v0, init, result_path)
save(result_optim, file = paste0(result_path, "result_optim.RData"))

num_cores <- n_cores
my_cluster <- parallel::makeCluster(num_cores, type = "PSOCK")
parallel::clusterEvalQ(my_cluster, {
  .libPaths("random_path")  # Adjust this to the correct path where R libraries are stored
})
doParallel::registerDoParallel(cl = my_cluster)
foreach::getDoParRegistered()

smoothed_theta_nomcmc <- vector("list", dim(result_optim[["mult_dir_samples"]])[3])

smoothed_theta_nomcmc <- foreach(c=1:dim(result_optim[["mult_dir_samples"]])[3], .packages = c("fenrir")) %dopar%
  {
    seed <- sample(1:10000, 1)
    res <- fenrir::fenrir_smooth(result_optim[["mult_dir_samples"]][,,c], F, G, gamma, W, M0, C0, Xi0, v0, observed_TT, N_total_list, seed)
    res$theta_smoothed
  }

smoothed_theta_nomcmc_matrix <- lapply(smoothed_theta_nomcmc, function(sample){t(do.call(rbind, sample))})

smoothed_theta_nomcmc_clr <- vector("list",length(smoothed_theta_nomcmc_matrix))

smoothed_theta_nomcmc_clr <- foreach(i=1:length(smoothed_theta_nomcmc), .packages = c("fido")) %dopar%
  {
    proportions <-  fido::alrInv(t(smoothed_theta_nomcmc_matrix[[i]]),1)
    return(t(fido::clr_array(proportions,2)))
  }

save(smoothed_theta_nomcmc, file = paste0(result_path,"result_theta_dir_without_mcmc.RData"))
save(smoothed_theta_nomcmc_clr, file = paste0(result_path,"result_theta_dir_without_mcmc_clr.RData"))

stopCluster(my_cluster)