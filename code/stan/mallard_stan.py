import os
import shutil
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

def clear_directory(directory):
    shutil.rmtree(directory)
    os.makedirs(directory)

def run_mcmc(D, Q, data_path, results_path, stan_model_path, infer_W = 0):
    # clear_directory(results_path)

    M0_list = [0,0,0,0]
    C0_list = [1,1,1,1]
    W_val = 0.5
    N_total_list = [693,693,693,693]
    N_obs_list = [158,158,158,158]
    Y = pd.read_csv(data_path+"Y_total.csv",header=None).to_numpy().astype(np.int64)
    observed_TT = pd.read_csv(data_path+"observed_indices.csv",header=None).to_numpy().flatten().astype(np.int64)

    q = Q 
    p = D-1  

    F = np.full((q, sum(N_total_list)), 1)
    G = [np.eye(q) * 1 for _ in range(sum(N_total_list))]
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
        "Y" : Y,
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
    fit = model.sample(data=data_dict,inits={"eta" : init}, output_dir=results_path, iter_sampling=100, iter_warmup=100, chains=2, parallel_chains=4)
    # Save the fit model using pickle
    with open(results_path + "fit.pkl", "wb") as f:
        pkl.dump(fit, f)


if __name__ == "__main__":

    cmdstanpy.set_cmdstan_path("/home/ayden/anaconda3/envs/stanpip/cmdstan-2.33.1")


    D = 10
    Q = 1
    data_path = "/home/ayden/Documents/Silverman Lab/code/fenrir_paper_code/data/mallard/"
    results_path = "/home/ayden/Documents/Silverman Lab/code/fenrir_paper_code/results/mallard/stan/"
    stan_model_dir = "/home/ayden/Documents/Silverman Lab/code/fenrir_paper_code/code/stan/"
    infer_W = 1

    stan_model_path = ""
    if infer_W:
        stan_model_path = stan_model_dir + "model_W.stan"
    else:
        stan_model_path = stan_model_dir + "model.stan"
    
    run_mcmc(D, Q, data_path, results_path, stan_model_path, infer_W)

    