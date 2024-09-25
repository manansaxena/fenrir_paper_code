import os
import cmdstanpy
from cmdstanpy import CmdStanModel
import numpy as np
import pandas as pd
import argparse
import pickle as pkl
import time
import csv

def run_mcmc(D, Q, data_path, results_path, stan_model_path, infer_W = 0):

    N_total_list = [673,673,673,673]
    N_obs_list = [135,133,133,136]
    Y_obs = pd.read_csv(data_path+"Y_obs.csv",header=None).to_numpy().astype(np.int64).transpose()
    observed_TT = pd.read_csv(data_path+"observed_indices.csv",header=None).to_numpy().flatten(order='F').astype(np.int64)

    q = Q 
    p = D-1  

    F = np.zeros((q, sum(N_total_list)))
    F[0, :] = 1
    G = [np.array([[1,1],[0,0.9]]) for _ in range(sum(N_total_list))]
    gamma = np.ones(sum(N_total_list))
    Xi0 = np.identity(p) * 10
    upsilon0 = D + 3
    init = np.zeros(shape=(p,sum(N_obs_list)))

    if not infer_W:
        W = []
    M0 = []
    C0 = []

    for i in range(len(N_total_list)):
        M0_val = np.zeros((q,p))
        C0_val = np.eye(q)
        M0.append(M0_val)
        C0.append(C0_val)
        if not infer_W:
            W.append([np.array([[0.3,0],[0,0.1]]) for _ in range(N_total_list[i])])

    if not infer_W:
        W = [item for sublist in W for item in sublist]

    if infer_W:
        inv_gamma_prior = {
            "a" : 30,
            "b" : 15,
            "c" : 30,
            "d" : 8
        }
    
    data_dict = {
        "N_total" : sum(N_total_list),
        "N_obs" : sum(N_obs_list),
        "num_timeseries": len(N_total_list),
        "observed_TT": observed_TT,
        "N_total_list": N_total_list,
        "p" : p,
        "q" : q,
        "Y_obs" : Y_obs,
        "FF" : F.T.tolist(),
        "GG" : G,
        "gamma" : gamma.tolist(),
        "M0" : M0,
        "C0" : C0,
        "Xi0" : Xi0,
        "upsilon0" : upsilon0,
    }

    if infer_W:
        data_dict = data_dict | inv_gamma_prior
    if not infer_W:
        data_dict = data_dict | {"WW" : W}

    # run mcmc
    model = CmdStanModel(stan_file=stan_model_path,cpp_options={'STAN_THREADS':'true'})
    start_time = time.time()
    fit = model.sample(data=data_dict,inits={"eta" : init}, output_dir=results_path, iter_sampling=3000, iter_warmup=1500, chains=4, parallel_chains=4) #3000,1500
    end_time = time.time()
    mcmc_time = end_time - start_time
    print("mcmc time", mcmc_time)

    # Save the fit model using pickle
    with open(results_path + "fit.pkl", "wb") as f:
        pkl.dump(fit, f)

def run_optimizer(D, Q, data_path, results_path, stan_model_path):

    N_total_list = [673,673,673,673]
    N_obs_list = [135,133,133,136]
    Y_obs = pd.read_csv(data_path+"Y_obs.csv",header=None).to_numpy().astype(np.int64).transpose()
    observed_TT = pd.read_csv(data_path+"observed_indices.csv",header=None).to_numpy().flatten(order='F').astype(np.int64)


    q = Q 
    p = D-1  

    F = np.zeros((q, sum(N_total_list)))
    F[0, :] = 1
    G = [np.array([[1,1],[0,0.9]]) for _ in range(sum(N_total_list))]
    gamma = np.ones(sum(N_total_list))
    Xi0 = np.identity(p) * 10
    upsilon0 = D + 3
    init = np.zeros(shape=(p,sum(N_obs_list)))

    W = []
    M0 = []
    C0 = []

    for i in range(len(N_total_list)):
        M0_val = np.zeros((q,p))
        C0_val = np.eye(q)
        M0.append(M0_val)
        C0.append(C0_val)
        W.append([np.array([[0.3,0],[0,0.1]]) for _ in range(N_total_list[i])])

    W = [item for sublist in W for item in sublist]
    
    data_dict = {
        "N_total" : sum(N_total_list),
        "N_obs" : sum(N_obs_list),
        "num_timeseries": len(N_total_list),
        "observed_TT": observed_TT,
        "N_total_list": N_total_list,
        "p" : p,
        "q" : q,
        "Y_obs" : Y_obs,
        "FF" : F.T.tolist(),
        "GG" : G,
        "WW" : W,
        "gamma" : gamma.tolist(),
        "M0" : M0,
        "C0" : C0,
        "Xi0" : Xi0,
        "upsilon0" : upsilon0,
    }

    iterations = [50,500,1000,10000]
    optim_time_list = []

    for iters in iterations:

        os.mkdir(results_path+str(iters))

        model = CmdStanModel(stan_file=stan_model_path)

        start_time = time.time()
        mle = model.optimize(data=data_dict,
                            seed=1,
                            inits={"eta" : init},
                            output_dir=results_path+str(iters),
                            algorithm="LBFGS",
                            iter = iters,
                            require_converged= True,
                            save_iterations=True,
                            refresh=1,
                            jacobian=True,
                            save_profile=True
                            )

        end_time = time.time()
        optimization_time = end_time - start_time

        optim_time_list.append([iters, optimization_time])

        for file in os.listdir(results_path+str(iters)):
            if file.endswith(".csv"):
                df = pd.read_csv(results_path+str(iters)+"/"+file, comment='#')
                header_index = df.columns.tolist().index('lp__')
                df = df.iloc[:, header_index:]
                df.to_csv(results_path+str(iters)+"/"+file, index=False)
    
    with open(results_path+'optimization_times.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Iteration", "OptimizationTime"])
        writer.writerows(optim_time_list)

if __name__ == "__main__":

    cmdstanpy.set_cmdstan_path("random_path" + ".conda/envs/stanpip/cmdstan-2.33.1") # set path where anaconda libraries are stored

    D = 10
    Q = 2
    data_path = "random_path_where_microbiome_data_are_stored"
    results_path = "random_path_where_results_on_microbiome_data_are_stored"
    stan_model_dir = "random_path_where_stan_model_are_stored"
    infer_W = 1

    stan_model_path = ""
    if infer_W:
        stan_model_path = stan_model_dir + "model_W_velocity.stan"
    else:
        stan_model_path = stan_model_dir + "model.stan"
    
    run_mcmc(D, Q, data_path, results_path, stan_model_path, infer_W)
    # run_optimizer(D, Q, data_path, results_path, stan_model_path)
    