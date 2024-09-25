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
Q <- 1
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
M0_list <- c(0,0,0,0)
C0_list <- c(1,1,1,1)

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


get_theta <- function(W_val, first_n_timeseries) {
  N_considered_list <- N_total_list[1:first_n_timeseries]
  N_considered <- sum(N_considered_list)
  observed_considered_TT <- observed_TT[1:N_considered]
  Y_considered_obs <- Y_obs[, 1:sum(observed_considered_TT)]

  W <- lapply(1:N_considered, function(x) diag(Q) * W_val)

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
  return(list(
    t(do.call(rbind, res$theta_smoothed)),
    t(do.call(rbind, res$theta0_smoothed)),
    res$Sigma
  ))
}

get_e <- function(theta, theta0, Sigma, N_total_list, first_n_timeseries) {
  N_considered_list <- N_total_list[1:first_n_timeseries]
  Sigma <- (Sigma + t(Sigma)) / 2.0
  e <- matrix(0, nrow = nrow(theta), ncol = ncol(theta))
  current_index <- 1
  for (i in seq_along(N_considered_list)) {
    N_current <- N_considered_list[i]
    for (j in 1:N_current) {
      index_global <- current_index + j - 1
      if (j == 1) {
        initial_theta <- theta0[,i]
        e[, index_global] <- theta[, index_global] - initial_theta
      } else {
        e[, index_global] <- theta[, index_global] - theta[, index_global - 1]
      }
    }
    current_index <- current_index + N_current
  }
  return(e)
}

get_W <- function(e, Sigma, nu_init, Xi_init, S) {
  Sigma <- (Sigma + t(Sigma)) / 2.0
  L <- t(chol((solve(Sigma))))
  e_transformed <- matrix(0, nrow = nrow(e), ncol = ncol(e))
  for (i in 1:nrow(e)) {
    e_transformed[i,] <- e[i,] %*% L
  }
  e_transformed = as.vector(e_transformed)
  nu_star <- nu_init + length(e_transformed)
  Xi_star <- (1/nu_star)*(nu_init * (Xi_init^2) + sum((e_transformed - mean(e_transformed))^2))
  W <- rinvgamma(S,nu_star/2,0.5 * nu_star * Xi_star)
  return(W)
}

time_list <- list()

z <- 4000 # number of iterations of gibbs sampler

theta_ct <- vector("list", num_timeseries)
e_ct <- vector("list", num_timeseries)
W_ct <- vector("list", num_timeseries)

num_cores <- parallel::detectCores() - 1
my_cluster <- parallel::makeCluster(num_cores, type = "PSOCK")
parallel::clusterEvalQ(my_cluster, {
  .libPaths("random_path")  # Adjust this to the correct path where R libraries are stored
})
doParallel::registerDoParallel(cl = my_cluster)
foreach::getDoParRegistered()

gibbs_start <- Sys.time()

# We run this gibbs sampler by considering different number of vessels.
results <- foreach(n_t = 1:num_timeseries, .packages = c("fenrir", "LaplacesDemon")) %dopar% {
  theta_z <- vector("list", z)
  e_z <- vector("list", z)
  W_z <- vector("list", z)
  for (i in 1:z) {
    if (i == 1) {
      W_z[[i]] <- rinvgamma(1,30,15)
      theta_z[[i]] <- get_theta(W_z[[i]], n_t)
    } else {
      theta_z[[i]] <- get_theta(W_z[[i - 1]], n_t)
    }
    e_z[[i]] <- get_e(theta_z[[i]][[1]], theta_z[[i]][[2]], theta_z[[i]][[3]][[1]], N_total_list, n_t)
    W_z[[i]] <- get_W(t(e_z[[i]]), theta_z[[i]][[3]][[1]], 60, 1/(2^0.5), 1)
  }

  list(theta_z = theta_z, e_z = e_z, W_z = W_z)
}

gibbs_end <- Sys.time()

for (n_t in 1:num_timeseries) {
  theta_ct[[n_t]] <- results[[n_t]]$theta_z
  e_ct[[n_t]] <- results[[n_t]]$e_z
  W_ct[[n_t]] <- results[[n_t]]$W_z
}

print("Gibbs sampling time: ")
gibbs_time <- gibbs_end - gibbs_start
print(gibbs_time)
time_list["gibbs_time"] <- gibbs_time

W_means <- numeric(num_timeseries)

for (i in 1:num_timeseries) {
  W_means[i] <- mean(unlist(lapply(W_ct[[i]][25:z], mean)))
}

save(W_ct, file = paste0(result_path, "W_gibbs.RData"))
save(W_means, file = paste0(result_path, "W_means_gibbs.RData"))

theta_smoothed <- array(0,dim=c(z,P,N_total))
for (i in 1:z){
  theta_smoothed[i,,] <- results[[1]][["theta_z"]][[i]][[1]]
}

smoothed_theta_nomcmc_clr <- array(0,dim=c(z,P+1,N_total))

for(c in 1:z){
  proportions <-  fido::alrInv(t(theta_smoothed[c,,]),1)
  smoothed_theta_nomcmc_clr[c,,] <- t(fido::clr_array(proportions,2))
}

save(theta_smoothed, file = paste0(result_path, "result_theta_dir_without_mcmc_gibbs.RData"))
save(smoothed_theta_nomcmc_clr, file = paste0(result_path, "result_theta_dir_without_mcmc_clr_gibbs.RData"))
save(time_list, file = paste0(result_path,"time_list.RData"))

stopCluster(my_cluster)