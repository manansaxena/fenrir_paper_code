import os
import cmdstanpy
from cmdstanpy import CmdStanModel
import numpy as np
import pandas as pd
import argparse
import pickle as pkl
import time
import csv

def replace_upper_tri(mat, val):
    mat[np.triu_indices_from(mat, k=1)] = val
    return mat

def run_mcmc(D, Q, data_path, results_path, stan_model_path, infer_W = 0):

    M0_list = pd.read_csv(data_path+"M0_list.csv",header=None).to_numpy().flatten().tolist()
    C0_list = pd.read_csv(data_path+"C0_list.csv",header=None).to_numpy().flatten().tolist()
    W_val = pd.read_csv(data_path+"W_val.csv",header=None).to_numpy().flatten().tolist()[0]
    N_total_list = pd.read_csv(data_path+"N_total_list.csv",header=None).to_numpy().flatten().tolist()
    N_obs_list = pd.read_csv(data_path+"N_obs_list.csv",header=None).to_numpy().flatten().tolist()
    Y_obs = pd.read_csv(data_path+"Y_obs.csv",header=None).to_numpy().astype(np.int64)
    observed_TT = pd.read_csv(data_path+"observed_indices.csv",header=None).to_numpy().flatten().astype(np.int64)

    q = Q 
    p = D-1  

    F = np.full((q, sum(N_total_list)), 1)
    G = [np.eye(q) * 1 for _ in range(sum(N_total_list))]
    gamma = np.ones(sum(N_total_list))
    Xi0 = np.identity(p) * 1
    upsilon0 = D + 3
    init = np.zeros(shape=(p,sum(N_obs_list)))

    if not infer_W:
        W = []
    M0 = []
    C0 = []

    for i in range(len(N_total_list)):
        M0_val = np.zeros((q,p))
        C0_val = np.eye(q)*C0_list[i]
        M0.append(replace_upper_tri(M0_val, M0_list[i]))
        C0.append(C0_val)
        if not infer_W:
            W.append([np.eye(q) * W_val for _ in range(N_total_list[i])])

    if not infer_W:
        W = [item for sublist in W for item in sublist]

    if infer_W:
        inv_gamma_prior = {
            "a" : 10,
            "b" : 5
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
    fit = model.sample(data=data_dict,inits={"eta" : init}, output_dir=results_path, iter_sampling=3000, iter_warmup=1500, chains=4, parallel_chains=4)
    end_time = time.time()
    mcmc_time = end_time - start_time
    print("mcmc time", mcmc_time)
    
    # Save the fit model using pickle
    with open(results_path + "fit.pkl", "wb") as f:
        pkl.dump(fit, f)

def optimize_gmdlm(D, Q, data_path, results_path, stan_model_path):

    M0_list = pd.read_csv(data_path+"M0_list.csv",header=None).to_numpy().flatten().tolist()
    C0_list = pd.read_csv(data_path+"C0_list.csv",header=None).to_numpy().flatten().tolist()
    W_val = pd.read_csv(data_path+"W_val.csv",header=None).to_numpy().flatten().tolist()[0]
    N_total_list = pd.read_csv(data_path+"N_total_list.csv",header=None).to_numpy().flatten().tolist()
    N_obs_list = pd.read_csv(data_path+"N_obs_list.csv",header=None).to_numpy().flatten().tolist()
    Y_obs = pd.read_csv(data_path+"Y_obs.csv",header=None).to_numpy().astype(np.int64)
    observed_TT = pd.read_csv(data_path+"observed_indices.csv",header=None).to_numpy().flatten().astype(np.int64)

    q = Q 
    p = D-1  

    F = np.full((q, sum(N_total_list)), 1)
    G = [np.eye(q) * 1 for _ in range(sum(N_total_list))]
    gamma = np.ones(sum(N_total_list))
    Xi0 = np.identity(p) * 1
    upsilon0 = D + 3
    init = np.zeros(shape=(p,sum(N_obs_list)))

    W = []
    M0 = []
    C0 = []

    for i in range(len(N_total_list)):
        M0_val = np.zeros((q,p))
        C0_val = np.eye(q)*C0_list[i]
        M0.append(replace_upper_tri(M0_val, M0_list[i]))
        C0.append(C0_val)
        W.append([np.eye(q) * W_val for _ in range(N_total_list[i])])

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

    parser = argparse.ArgumentParser(description='for stan collapsed model')
    parser.add_argument('D', type=int)
    parser.add_argument('Q', type=int)
    parser.add_argument('data_path', type=str)
    parser.add_argument('results_path', type=str)
    parser.add_argument('stan_model_dir', type=str)
    parser.add_argument('infer_W', type=int)
    parser.add_argument('run_mcmc',type=int)

    args = parser.parse_args()

    stan_model_path = ""
    if args.infer_W:
        stan_model_path = args.stan_model_dir + "model_W.stan"
    else:
        stan_model_path = args.stan_model_dir + "model.stan"
    
    if args.run_mcmc:
        run_mcmc(args.D, args.Q, args.data_path, args.results_path, stan_model_path, args.infer_W)
    else:
        optimize_gmdlm(args.D, args.Q, args.data_path, args.results_path, stan_model_path)

    