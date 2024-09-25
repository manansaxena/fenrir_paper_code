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

F <- matrix(0, nrow = Q, ncol = N_total)
F[1, ] <- 1

G <- lapply(1:N_total, function(x) matrix(c(1,0,1,0.9),nrow=Q,ncol=Q))
gamma <- rep(1, N_total)
M0 <- lapply(1:num_timeseries, function(x) matrix(0, nrow = Q, ncol = P))
C0 <- lapply(1:num_timeseries, function(x) diag(Q))
Xi0 <- 10 * diag(P)
v0 <- D + 3
init <- matrix(0, nrow = P, ncol = N_obs)


get_theta <- function(W_theta_val, W_alpha_val, first_n_timeseries) {
  N_considered_list <- N_total_list[1:first_n_timeseries]
  N_considered <- sum(N_considered_list)
  observed_considered_TT <- observed_TT[1:N_considered]
  Y_considered_obs <- Y_obs[, 1:sum(observed_considered_TT)]

  W_val <- matrix(c(W_theta_val,0,0,W_alpha_val),nrow=Q,ncol=Q)
  W <- lapply(1:N_considered, function(x) W_val)

  result_optim <- optimizer(
    Y_considered_obs, observed_considered_TT, N_considered_list,
    F[, 1:N_considered, drop = FALSE], G[1:N_considered], gamma[1:N_considered], W,
    M0[1:first_n_timeseries], C0[1:first_n_timeseries], Xi0, v0, init[, 1:sum(observed_considered_TT)],
    result_path
  )
  seed <- sample(1:10000, 1)
  which_multdir <- sample(1:dim(result_optim[["mult_dir_samples"]])[3],1)
  res <- fenrir::fenrir_smooth(
      result_optim[["mult_dir_samples"]][,,which_multdir], F[, 1:N_considered, drop = FALSE], G[1:N_considered], gamma[1:N_considered], W,
      M0[1:first_n_timeseries], C0[1:first_n_timeseries], Xi0, v0,
      observed_considered_TT, N_considered_list, seed
  )
  combined_matrix <- matrix(nrow = P, ncol = N_considered*Q)
  column_index <- 1
  for (j in 1:N_considered) {
    current_matrix <- res[["theta_smoothed"]][[j]]
    combined_matrix[, column_index] <- t(current_matrix)[,1]
    combined_matrix[, N_considered+column_index] <- t(current_matrix)[,2]
    column_index <- column_index + 1
  }

  return(list(
    combined_matrix,
    res$theta0_smoothed,
    res$Sigma
  ))
}

get_e <- function(theta, theta0, N_total_list, first_n_timeseries) {
  N_considered_list <- N_total_list[1:first_n_timeseries]
  N_considered <- sum(N_considered_list)
  e_theta <- matrix(0, nrow = D-1, ncol = N_considered)
  e_alpha <- matrix(0, nrow = D-1, ncol = N_considered)

  current_index <- 1
  for (i in seq_along(N_considered_list)) {
    N_current <- N_considered_list[i]
    for (j in 1:N_current) {
      index_global <- current_index + j - 1
      if (j == 1) {
        initial_theta <- theta0[[i]][1,]
        initial_alpha <- theta0[[i]][2,]
        e_theta[, index_global] <- theta[, index_global] - initial_theta - initial_alpha
        e_alpha[, index_global] <- theta[, index_global + N_considered] - 0.9 * initial_alpha
      } else {
        e_theta[, index_global] <- theta[, index_global] - theta[, index_global - 1] - theta[, index_global + N_considered - 1]
        e_alpha[, index_global] <- theta[, index_global + N_considered] - 0.9 * theta[, index_global + N_considered - 1]
      }
    }
    current_index <- current_index + N_current
  }
  return(list(e_theta, e_alpha))
}

get_W <- function(e_theta, e_alpha, Sigma, nu_init_theta, Xi_init_theta, nu_init_alpha, Xi_init_alpha, S) {
  Sigma <- (Sigma + t(Sigma)) / 2.0
  L <- t(chol((solve(Sigma))))
  e_theta_transformed <- matrix(0, nrow = nrow(e_theta), ncol = ncol(e_theta))
  e_alpha_transformed <- matrix(0, nrow = nrow(e_alpha), ncol = ncol(e_alpha))
  for (i in 1:nrow(e_theta)) {
    e_theta_transformed[i,] <- e_theta[i,] %*% L
    e_alpha_transformed[i,] <- e_alpha[i,] %*% L
  }
  e_theta_transformed = as.vector(e_theta_transformed)
  e_alpha_transformed = as.vector(e_alpha_transformed)

  nu_star_theta <- nu_init_theta + length(e_theta_transformed)
  Xi_star_theta <- (1/nu_star_theta)*(nu_init_theta * (Xi_init_theta^2) + sum((e_theta_transformed - mean(e_theta_transformed))^2))
  W_theta <- rinvgamma(S,nu_star_theta/2,0.5 * nu_star_theta * Xi_star_theta)

  nu_star_alpha <- nu_init_alpha + length(e_alpha_transformed)
  Xi_star_alpha <- (1/nu_star_alpha)*(nu_init_alpha * (Xi_init_alpha^2) + sum((e_alpha_transformed - mean(e_alpha_transformed))^2))
  W_alpha <- rinvgamma(S,nu_star_alpha/2,0.5 * nu_star_alpha * Xi_star_alpha)
  return(list(W_theta, W_alpha))
}

time_list <- list()

z <- 4000 # number of iterations of gibbs sampler

gibbs_start <- Sys.time()

theta_z <- vector("list", z)
e_z <- vector("list", z)
W_z <- vector("list", z)
n_t <- 4

# Unlike mallard_gibbs.R, I run this model considering all four vessels.
for (i in 1:z) {
  W_z[[i]] <- numeric(2)

  if (i == 1) {
    W_z[[i]][1] <- rinvgamma(1, 30, 15)
    W_z[[i]][2] <- rinvgamma(1, 30, 8)
    theta_z[[i]] <- get_theta(W_z[[i]][1], W_z[[i]][2], n_t)
  } else {
    theta_z[[i]] <- get_theta(W_z[[i - 1]][1], W_z[[i - 1]][2], n_t)
  }

  e_z[[i]] <- get_e(theta_z[[i]][[1]], theta_z[[i]][[2]], N_total_list, n_t)
  W_comb <- get_W(t(e_z[[i]][[1]]), t(e_z[[i]][[2]]), theta_z[[i]][[3]][[1]], 60, 1/(2^0.5), 60, (4/15)^0.5, 1)
  W_z[[i]][1] <- W_comb[[1]]
  W_z[[i]][2] <- W_comb[[2]]
}

gibbs_end <- Sys.time()

print("Gibbs sampling time: ")
gibbs_time <- gibbs_end - gibbs_start
print(gibbs_time)
time_list["gibbs_time"] <- gibbs_time

W_means <- numeric(2)
for (val in W_z) {
  W_means <- W_means + val
}
W_means <- W_means / length(W_z)

save(W_z, file = paste0(result_path, "W_gibbs.RData"))
save(W_means, file = paste0(result_path, "W_means_gibbs.RData"))

theta_smoothed <- array(0,dim=c(z,P,N_total*Q))
for (i in 1:z){
  theta_smoothed[i,,] <- theta_z[[i]][[1]]
}

smoothed_theta_nomcmc_clr <- array(0,dim=c(z,P+1,N_total*Q))

for(c in 1:z){
  proportions <-  fido::alrInv(t(theta_smoothed[c,,]),1)
  smoothed_theta_nomcmc_clr[c,,] <- t(fido::clr_array(proportions,2))
}

save(theta_smoothed, file = paste0(result_path, "result_theta_dir_without_mcmc_gibbs.RData"))
save(smoothed_theta_nomcmc_clr, file = paste0(result_path, "result_theta_dir_without_mcmc_clr_gibbs.RData"))
save(time_list, file = paste0(result_path,"time_list.RData"))