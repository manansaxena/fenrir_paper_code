# set the path where all the R libraries are stored.
.libPaths(c("random_path", .libPaths()))

library(fenrir)
library(fido)
library(coda)
library(abind)
library(driver)
library(LaplacesDemon)
library(compositions)
library(doParallel)
library(foreach)

optimizer <- function(Y_obs, observed_TT, N_total_list, F, G, gamma, W, M0, C0, Xi0, v0, init, result_path) {
  start_time <- Sys.time()
  res <- fenrir::fenrir_optim(
    Y_obs = Y_obs, observed_TT = observed_TT, N_total_list = N_total_list, F = F, G = G, gamma = gamma, W = W, M0 = M0, C0 = C0, Xi0 = Xi0, v0 = v0, init = init,
    log_probs_path = paste0(result_path, "logprobs_cpp.csv"), num_dirsamples = 2000, pseudocount = 0.5, max_iter = 10000, eps_f = 1e-7, eps_g = 1e-6
  )
  end_time <- Sys.time()
  res$total_time <- end_time - start_time
  return(res)
}

D <- 10
Q <- 2
data_path <- "random_path_where_microbiome_data_are_stored" # specify the path where microbiome(mallard) data is stored
result_path <- "random_path_where_results_on_microbiome_data_are_stored" # specify the path where results would be stored

# Load data
data <- readRDS(paste0(data_path, "data.rds"))

Y_obs <- as.matrix(data$Y_obs)
rownames(Y_obs) <- NULL
colnames(Y_obs) <- NULL
Y_obs <- t(Y_obs)
observed_TT <- as.numeric(as.vector(data$observed_TT))

N_total_list <- c(673,673,673,673)
N_total <- length(observed_TT)
N_obs <- dim(Y_obs)[2]
P <- D - 1
num_timeseries <- length(N_total_list)
W_val <- matrix(c(0.3,0,0,0.1),nrow=Q,ncol=Q)
W <- lapply(1:N_total, function(x) W_val)
F <- matrix(0, nrow = Q, ncol = N_total)
F[1, ] <- 1
G <- lapply(1:N_total, function(x) matrix(c(1,0,1,0.9),nrow=Q,ncol=Q))
gamma <- rep(1, N_total)
M0 <- lapply(1:num_timeseries, function(x) matrix(0, nrow = Q, ncol = P))
C0 <- lapply(1:num_timeseries, function(x) diag(Q))
Xi0 <- 10 * diag(P)
v0 <- D + 3
init <- matrix(0, nrow = P, ncol = N_obs)

result_optim <- optimizer(Y_obs, observed_TT, N_total_list, F, G, gamma, W, M0, C0, Xi0, v0, init, result_path)
save(result_optim, file = paste0(result_path, "result_optim.RData"))

smoothed_theta_nomcmc <- vector("list", dim(result_optim[["mult_dir_samples"]])[3])

for (c in 1:dim(result_optim[["mult_dir_samples"]])[3]) {
  seed <- sample(1:10000, 1)
  res <- fenrir::fenrir_smooth(result_optim[["mult_dir_samples"]][,,c], F, G, gamma, W, M0, C0, Xi0, v0, observed_TT, N_total_list, seed)
  smoothed_theta_nomcmc[[c]] <- res$theta_smoothed
}

smoothed_theta_nomcmc_matrix <- array(0, dim=c(length(smoothed_theta_nomcmc), Q, P, N_total))
for (i in 1:length(smoothed_theta_nomcmc)) {
  for (j in 1:N_total){
    smoothed_theta_nomcmc_matrix[i,,,j] = smoothed_theta_nomcmc[[i]][[j]]
  }
}

smoothed_theta_nomcmc_clr <- array(0, dim=c(length(smoothed_theta_nomcmc),Q,P+1,N_total))

for(i in 1:length(smoothed_theta_nomcmc)){
  for (j in 1:Q){
    proportions <-  fido::alrInv(t(smoothed_theta_nomcmc_matrix[i,j,,]),1)
    smoothed_theta_nomcmc_clr[i,j,,] <- t(fido::clr_array(proportions,2))
  }
}

save(smoothed_theta_nomcmc, file = paste0(result_path, "result_theta_dir_without_mcmc.RData"))
save(smoothed_theta_nomcmc_clr, file = paste0(result_path, "result_theta_dir_without_mcmc_clr.RData"))
