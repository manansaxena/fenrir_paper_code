library(fenrir)
library(fido)
library(coda)
library(abind)
library(driver)
library(compositions)

optimizer <- function(Y_obs, observed_TT, N_total_list, F, G, gamma, W, M0, C0, Xi0, v0, init, result_path) {
  start_time <- Sys.time()
  res <- fenrir::fenrir_optim(
    Y_obs = Y_obs, observed_TT = observed_TT, N_total_list = N_total_list, F = F, G = G, gamma = gamma, W = W, M0 = M0, C0 = C0, Xi0 = Xi0, v0 = v0, init = init,
    log_probs_path = paste0(result_path, "logprobs_cpp.csv"), num_dirsamples = 10, pseudocount = 0.5, max_iter = 10000, eps_f = 1e-6, eps_g = 1e-5
  )
  end_time <- Sys.time()
  res$total_time <- end_time - start_time
  return(res)
}

D <- 10
Q <- 1
data_path <- "/home/ayden/Documents/Silverman Lab/code/fenrir_paper_code/data/mallard/"
result_path <- "/home/ayden/Documents/Silverman Lab/code/fenrir_paper_code/results/mallard/cpp/"
run_mcmc <- 1
n_cores <- 2

# Load data
data <- readRDS(paste0(data_path, "data.rds"))

Y_obs <- data$Y_obs_combined
Y_obs <- aperm(Y_obs, c(3, 1, 2))
Y_obs <- matrix(Y_obs, nrow = 10, ncol = 158 * 4)
observed_TT <- as.numeric(as.vector(data$observed_indices_combined))

N_total_list <- c(693,693,693,693)
N_total <- length(observed_TT)
N_obs <- dim(Y_obs)[2]
P <- D - 1
num_timeseries <- length(N_total_list)
M0_list <- c(0,0,0,0)
C0_list <- c(1,1,1,1)
W_val <- 0.5
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
Xi0 <- 10 * diag(P)
v0 <- D + 3
init <- matrix(0, nrow = P, ncol = N_obs)

result_optim <- optimizer(Y_obs, observed_TT, N_total_list, F, G, gamma, W, M0, C0, Xi0, v0, init, result_path)
saveRDS(result_optim, file = paste0(result_path, "result_optim.rds"))


smoothed_theta_nomcmc <- vector("list", dim(result_optim[["mult_dir_samples"]])[3])

for (c in 1:dim(result_optim[["mult_dir_samples"]])[3]) {
  seed <- sample(1:10000, 1)
  res <- fenrir::fenrir_smooth(result_optim[["mult_dir_samples"]][,,c], F, G, gamma, W, M0, C0, Xi0, v0, observed_TT, N_total_list, seed)
  smoothed_theta_nomcmc[[c]] <- res$theta_smoothed
}

smoothed_theta_nomcmc_matrix <- lapply(smoothed_theta_nomcmc, function(sample){t(do.call(rbind, sample))})

smoothed_theta_nomcmc_clr <- vector("list",length(smoothed_theta_nomcmc_matrix))
for(i in 1:length(smoothed_theta_nomcmc)) {
  proportions <-  fido::alrInv(t(smoothed_theta_nomcmc_matrix[[i]]),1)
  smoothed_theta_nomcmc_clr[[i]] <- t(fido::clr_array(proportions,2))
}

saveRDS(smoothed_theta_nomcmc, file = paste0(result_path,"result_theta_dir_without_mcmc.rds"))
saveRDS(smoothed_theta_nomcmc_clr, file = paste0(result_path,"result_theta_dir_without_mcmc_clr.rds"))


if(run_mcmc == 1){
  etadir_matrix <- apply(result_optim[["mult_dir_samples"]], 3, function(x) as.vector(t(x)))
  etadir_chol <- chol(cov(t(etadir_matrix)))
  num_chains <- 10
  num_steps <- 100
  seed <- sample(1:10000, 1)
  scale_factor_of_proposal <- 0.008

  selected_indices <- sample(dim(result_optim[["mult_dir_samples"]])[3], num_chains)
  initial_states <- matrix(nrow = sum(N_obs) * P, ncol = num_chains)

  for (i in 1:num_chains) {
    initial_states[, i] <- as.vector(t(result_optim[["mult_dir_samples"]][, , selected_indices[i]]))
  }

  mcmc_res <- fenrir::fenrir_mh(
    Y_obs = Y_obs, observed_TT = observed_TT, N_total_list = N_total_list, N_obs = N_obs, F = F, G = G, gamma = gamma, W = W, M0 = M0, C0 = C0, Xi0 = Xi0, v0 = v0, 
    initial_states = initial_states, run_mcmc = 1, num_chains = num_chains, num_steps = num_steps, seed = seed, n_cores = n_cores, scale_factor_of_proposal = scale_factor_of_proposal,
    cholesky = etadir_chol
  )

  saveRDS(mcmc_res, file = paste0(result_path, "mcmc_res.rds"))

  smoothed_theta <- list()

  for (s in 1:length(mcmc_res[["particles"]])) {
    param_samples_matrix <- mcmc_res[["particles"]][[s]]
    smoothed_theta[[s]] <- list()
    for (i in 1:ncol(param_samples_matrix)) {
      param_samples <- t(matrix(param_samples_matrix[, i], nrow = N_obs, ncol = P))
      seed <- sample(1:10000, 1)
      res <- fenrir::fenrir_smooth(param_samples, F, G, gamma, W, M0, C0, Xi0, v0, observed_TT, N_total_list, seed)
      smoothed_theta[[s]][[i]] <- res$theta_smoothed
    }
  }
  
  smoothed_theta_clr <- list()
  for (s in 1:length(smoothed_theta)){
    smoothed_theta_matrix <- lapply(smoothed_theta[[s]], function(sample){t(do.call(rbind, sample))})
    smoothed_theta_clr_sample <- vector("list",length(smoothed_theta_matrix))
    for(i in 1:length(smoothed_theta[[s]])) {
      proportions <-  fido::alrInv(t(smoothed_theta_matrix[[i]]),1)
      smoothed_theta_clr_sample[[i]] <- t(fido::clr_array(proportions,2))
    }
    smoothed_theta_clr[[s]] <- smoothed_theta_clr_sample
  }
  
  saveRDS(smoothed_theta, file = paste0(result_path, "result_theta_dir_with_mcmc.rds"))
  saveRDS(smoothed_theta_clr, file = paste0(result_path, "result_theta_dir_with_mcmc_clr.rds"))
}
