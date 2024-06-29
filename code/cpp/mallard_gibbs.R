library(fenrir)
library(fido)
library(coda)
library(abind)
library(driver)
library(compositions)
library(CholWishart)
library(LaplacesDemon)
library(doParallel)
library(foreach)

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
run_mcmc <- 0
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


get_theta <- function(W_val, first_n_timeseries, one_theta = TRUE) {
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

  if (one_theta) {
    seed <- sample(1:10000, 1)
    res <- fenrir::fenrir_smooth(
      result_optim[["optim_eta"]], F[, 1:N_considered, drop = FALSE], G[1:N_considered], gamma[1:N_considered], W,
      M0[1:first_n_timeseries], C0[1:first_n_timeseries], Xi0, v0,
      observed_considered_TT, N_considered_list, seed
    )
    return(list(
      t(do.call(rbind, res$theta_smoothed)),
      t(do.call(rbind, res$theta0_smoothed)),
      res$Sigma
    ))
  } else {
    smoothed_theta_nomcmc <- vector("list", dim(result_optim[["mult_dir_samples"]])[3])
    for (c in 1:dim(result_optim[["mult_dir_samples"]])[3]) {
      seed <- sample(1:10000, 1)
      res <- fenrir::fenrir_smooth(
        result_optim[["mult_dir_samples"]][, , c], F[, 1:N_considered, drop = FALSE], G[1:N_considered], gamma[1:N_considered], W,
        M0[1:first_n_timeseries], C0[1:first_n_timeseries], Xi0, v0,
        observed_considered_TT, N_considered_list, seed
      )

      smoothed_theta_nomcmc[[c]] <- res$theta_smoothed
    }
    return(smoothed_theta_nomcmc)
  }
}

# theta : D-1 x N, M0: 1 X D-1 list, C0: 1 x 1 list, Sigma: D-1xD-1
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
  L <- t(chol(solve(Sigma)))
  e_transformed <- as.vector(e %*% t(L))
  nu_star <- nu_init + length(e_transformed)
  Xi_star <- Xi_init + (t(e_transformed) - c(mean(e_transformed))) %*% t((t(e_transformed) - c(mean(e_transformed))))
  # W <- rInvWishart(S, nu_star, Xi_star)[,,]
  W <- rinvgamma(S,nu_star,Xi_star)
  return(W)
}


z <- 200 # number of iterations of gibbs sampler

theta_ct <- vector("list", num_timeseries)
e_ct <- vector("list", num_timeseries)
W_ct <- vector("list", num_timeseries)

# num_cores <- parallel::detectCores() - 1
num_cores <- n_cores
my_cluster <- parallel::makeCluster(num_cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my_cluster)
foreach::getDoParRegistered()

gibbs_start <- Sys.time()

results <- foreach(n_t = 1:num_timeseries, .packages = c("fenrir", "CholWishart", "LaplacesDemon")) %dopar% {
  theta_z <- vector("list", z)
  e_z <- vector("list", z)
  W_z <- vector("list", z)

  for (i in 1:z) {
    if (i == 1) {
      # W_z[[i]] <- rInvWishart(1, 5, 0.5)[1]
      W_z[[i]] <- rinvgamma(1,10,5)
      theta_z[[i]] <- get_theta(W_z[[i]], n_t, TRUE)
    } else {
      theta_z[[i]] <- get_theta(mean(W_z[[i - 1]]), n_t, TRUE)
    }
    e_z[[i]] <- get_e(theta_z[[i]][[1]], theta_z[[i]][[2]], theta_z[[i]][[3]][[1]], N_total_list, n_t)
    W_z[[i]] <- get_W(t(e_z[[i]]), theta_z[[i]][[3]][[1]], 10, 5, 1000)
  }

  list(theta_z = theta_z, e_z = e_z, W_z = W_z)
}

gibbs_end <- Sys.time()

stopCluster(my_cluster)

# Organize results back into the lists
for (n_t in 1:num_timeseries) {
  theta_ct[[n_t]] <- results[[n_t]]$theta_z
  e_ct[[n_t]] <- results[[n_t]]$e_z
  W_ct[[n_t]] <- results[[n_t]]$W_z
}

# End timing
print("Gibbs sampling time: ")
print(gibbs_end - gibbs_start)

W_means <- numeric(num_timeseries)

# Compute the mean for each element in the list of lists
for (i in 1:num_timeseries) {
  W_means[i] <- mean(unlist(lapply(W_ct[[i]][(z/2):z], mean)))
}

smoothed_theta_nomcmc <- get_theta(W_means[length(W_means)], num_timeseries, FALSE)

smoothed_theta_nomcmc_matrix <- lapply(smoothed_theta_nomcmc, function(sample) {
  t(do.call(rbind, sample))
})

smoothed_theta_nomcmc_clr <- vector("list", length(smoothed_theta_nomcmc_matrix))
for (i in 1:length(smoothed_theta_nomcmc)) {
  proportions <- fido::alrInv(t(smoothed_theta_nomcmc_matrix[[i]]), 1)
  smoothed_theta_nomcmc_clr[[i]] <- t(fido::clr_array(proportions, 2))
}

saveRDS(W_means, file = paste0(result_path, "W_means_gibbs.rds"))
saveRDS(smoothed_theta_nomcmc, file = paste0(result_path, "result_theta_dir_without_mcmc_gibbs.rds"))
saveRDS(smoothed_theta_nomcmc_clr, file = paste0(result_path, "result_theta_dir_without_mcmc_clr_gibbs.rds"))
