# set the path where all the R libraries are stored.
.libPaths(c("random_path", .libPaths()))

library(fido)
library(MCMCpack)

simulate_data_for_multiple_ts <- function(D, Q, rseed, W_val, percent_of_missing, num_timeseries, N_total_list, M0_list, C0_list){
  Q <- Q
  D <- D
  P <- D-1
  # we have one upsilon and Xi for the entire long time series.
  upsilon <- D + 3
  Xi0 <- diag(P)
  Sigma <- riwish(upsilon, Xi0)

  Theta_total <- vector("list", num_timeseries)
  Y_total <- vector("list", num_timeseries)
  Y_obs_total <- vector("list", num_timeseries)
  observed_indices_total <- vector("list", num_timeseries)
  N_obs_list <- vector("list", num_timeseries)

  # The rest varies for each time series. M and C are reinitialized for each time series.
  for(timeseries in c(1:num_timeseries)){
    num_missing_values <- floor(0.01 * percent_of_missing * N_total_list[timeseries])
    selected_values <- sample(1:N_total_list[timeseries], num_missing_values)

    N <- N_total_list[timeseries] - num_missing_values
    N_obs_list[timeseries] <- N
    M0 <- matrix(0, nrow = Q, ncol = P)
    M0[upper.tri(M0)] <- M0_list[timeseries]
    C0 <- diag(Q) * C0_list[timeseries]

    Theta0 <- M0 + t(chol(C0))%*%matrix(rnorm(Q*(P)), nrow=Q)%*%chol(Sigma)

    gamma <- rep(1,N)
    W_val <- W_val
    W <- lapply(1:N, function(x) diag(Q) * W_val)
    W <- lapply(W, function(x) matrix(x, nrow = Q, ncol = Q))
    G <- lapply(1:N, function(x) diag(Q) * 1)
    G <- lapply(G, function(x) matrix(x, nrow = Q, ncol = Q))
    F <- matrix(1, nrow = Q, ncol = N)
    Theta <- lapply(1:N, function(x) matrix(0, nrow = Q, ncol = P))
    Theta <- lapply(Theta, function(x) matrix(x, nrow = Q, ncol = P))
    eta <- matrix(0, nrow = P, ncol = N)

    for(i in c(1:N)){
      if(i == 1){
        ohm <- matrix(0, nrow = Q, ncol = P) + t(chol(W[[1]]))%*%matrix(rnorm(Q*(P)), nrow=Q)%*%chol(Sigma)
        v <- numeric(P) + t(chol(gamma[i]*Sigma))%*%matrix(rnorm(P), nrow = P)
        Theta[[1]] <- G[[1]] %*% Theta0 + ohm
        eta[,1] <- t(t(F[,1]) %*% Theta[[1]] + t(v))
      }else{
        ohm <- matrix(0, nrow = Q, ncol = P) + t(chol(W[[i]]))%*%matrix(rnorm(Q*(P)), nrow=Q)%*%chol(Sigma)
        v <- numeric(P) + t(chol(gamma[i]*Sigma))%*%matrix(rnorm(P), nrow = P)
        Theta[[i]] <- G[[i]] %*% Theta0 + ohm
        eta[,i] <- t(t(F[,i]) %*% Theta[[i]] + t(v))
      }
    }
    Pi <- t(alrInv(t(eta)))
    Y <- matrix(0, D, N)
    for (i in 1:N) Y[,i] <- rmultinom(1, sample(0:5000), prob = Pi[,i])
    Y_missing <- matrix(-1, D, N_total_list[timeseries])
    observed_indices <- rep(0,N_total_list[timeseries])
    j <- 1
    for (i in 1:N_total_list[timeseries]) {
      if (!(i %in% selected_values)) {
        Y_missing[, i] <- Y[, j]
        observed_indices[i] <- 1
        j <- j + 1
      }
    }
    Theta_total[[timeseries]] <- Theta
    Y_total[[timeseries]] <- Y_missing
    observed_indices_total[[timeseries]] <- observed_indices
    Y_obs_total[[timeseries]] <- Y
  }
  return (list(Theta_total = Theta_total, Y_total = Y_total, observed_indices_total = observed_indices_total, Y_obs_total = Y_obs_total, N_obs_list = N_obs_list))
}

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 9){
  stop(paste("Usage: Rscript generate_simulated_data.R {D} {Q} {rseed} {W_val} {percent_of_missing} {N_total_list} {M0_list} {C0_list} {out_dir}"))
}

D <- as.integer(args[1])
Q <- as.integer(args[2])
rseed <- as.integer(args[3])
W_val <- as.numeric(args[4])
percent_of_missing <- as.numeric(args[5])
N_total_list <- as.integer(unlist(strsplit(args[6], ",")))
M0_list <- as.numeric(unlist(strsplit(args[7], ",")))
C0_list <- as.numeric(unlist(strsplit(args[8], ",")))
out_dir <- paste0(args[9], "/N", sum(N_total_list), "_D", D, "_Q", Q, "_R", rseed, "_W", W_val, "_pm", percent_of_missing, "/")

set.seed(rseed)
num_timeseries <- length(N_total_list)

data <- simulate_data_for_multiple_ts(D, Q, rseed, W_val, percent_of_missing, num_timeseries, N_total_list, M0_list, C0_list)

Y_combined <- do.call(cbind, data[["Y_total"]])
observed_indices_combined <- unlist(data[["observed_indices_total"]])
Y_obs_combined <- do.call(cbind, data[["Y_obs_total"]])
Theta_flat <- unlist(data[["Theta_total"]])
Theta_matrix <- matrix(Theta_flat,nrow=dim(Y_obs_combined)[1]-1,ncol=dim(Y_obs_combined)[2])
N_obs_list <- do.call(cbind, data[["N_obs_list"]])

dir.create(out_dir, showWarnings = FALSE)
write.table(Y_combined, paste0(out_dir,"Y.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(observed_indices_combined, paste0(out_dir,"observed_indices.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(Y_obs_combined, paste0(out_dir,"Y_obs.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(Theta_matrix, paste0(out_dir,"Theta.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(N_total_list, paste0(out_dir,"N_total_list.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(M0_list, paste0(out_dir,"M0_list.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(C0_list, paste0(out_dir,"C0_list.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(W_val, paste0(out_dir,"W_val.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(percent_of_missing, paste0(out_dir,"percent_of_missing.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(N_obs_list, paste0(out_dir,"N_obs_list.csv"), row.names = FALSE, col.names = FALSE, sep = ",")


data_saved <- list(Y_combined = Y_combined,
                   observed_indices_combined = observed_indices_combined,
                   Y_obs_combined = Y_obs_combined,
                   Theta_matrix = Theta_matrix,
                   N_total_list = N_total_list,
                   M0_list = M0_list,
                   C0_list = C0_list,
                   W_val = W_val,
                   percent_of_missing = percent_of_missing,
                   N_obs_list = N_obs_list)
saveRDS(data_saved, file=paste0(out_dir, "data.rds"))

