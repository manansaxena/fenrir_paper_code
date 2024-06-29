import sys
import os
import subprocess
import time
import shutil
import numpy as np

def clear_directory(directory):
    shutil.rmtree(directory)
    os.makedirs(directory)

def generate_params(N_vals):
    # each time series would be of 100 timepoints
    N_total_list = []
    M0_list = []
    C0_list = []
    num_timeseries = []
    for N_val in N_vals:
        num_timeseries.append(int(int(N_val)/100))
        N_total_list.append(np.repeat(100, int(int(N_val)/100)))
        M0_list.append(np.random.uniform(0.1, 1.0, int(int(N_val)/100)))
        C0_list.append(np.random.uniform(1, 1.5, int(int(N_val)/100)))
    return N_total_list, M0_list, C0_list, num_timeseries

if len(sys.argv) < 2:
    print("Usage: python generate_runtime_jobs.py {vary} {cores}")
    sys.exit()

vary = sys.argv[1]
cores = sys.argv[2]
model_dir = 'batch_jobs'

if vary == 'N':
    N_vals = ['200','300']
    D_vals = ['3']
    Q_vals = ['1']
    
elif vary == 'D':
    N_vals = ['500']
    D_vals = ['3', '25', '50', '75']
    Q_vals = ['1']

N_total_list, M0_list, C0_list, num_timeseries = generate_params(N_vals)
W_val = "0.45"
percent_of_missing = "5"
infer_W = "1"
run_mcmc = "1"

methods = ['mycpp']
rseed = ['1','2']

print(f"Generating {len(N_vals)*len(D_vals)*len(Q_vals)*len(methods)*len(rseed)} slurm scripts...")

simulation_ret = -1
sys_response = ''

file_suffix = ""
file_suffix = file_suffix + "_" + cores

for N in range(0,len(N_vals)):
    for D in range(0,len(D_vals)):
        for Q in range(0,len(Q_vals)):
            for seed in range(0,len(rseed)):

                N_total_str = ','.join(map(str, N_total_list[N]))
                M0_str = ','.join(map(str, M0_list[N]))
                C0_str = ','.join(map(str, C0_list[N]))
                simulation_ret = subprocess.call(["Rscript", "generate_simulated_data.R", D_vals[D], Q_vals[Q], rseed[seed], W_val, percent_of_missing, N_total_str, M0_str, C0_str, "../data/simulated"])
                if simulation_ret != 0:
                    print(f"Problem generating simulated data for {N_vals[N]}, {D_vals[D]}, {Q_vals[Q]} case (varying {vary})!")
                    exit(1)

                for m_idx in methods:

                    filename = f'../scripts/{vary}-varying_N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}_{m_idx}_{infer_W}_{run_mcmc}.sh'
                    
                    with open(filename, 'w') as fh:
                        if not os.path.exists(f'../results/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}_{m_idx}_{infer_W}_{run_mcmc}'):
                                os.mkdir(f'../results/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}_{m_idx}_{infer_W}_{run_mcmc}')

                        if m_idx == "sc":
                            fh.write('#!/bin/bash\n')
                            fh.write(f'cd "/home/ayden/Documents/Silverman Lab/code/fenrir_paper_code/" \n\n')
                            fh.write(f'/home/ayden/anaconda3/envs/stanpip/bin/python code/stan/run_stan.py {D_vals[D]} {Q_vals[Q]} "data/simulated/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}/" "results/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}_{m_idx}_{infer_W}_{run_mcmc}/" "code/stan/" {infer_W} {run_mcmc}')

                        elif m_idx == "mycpp": 

                            result = subprocess.run(['Rscript', '-e', 'if (!"fenrir" %in% installed.packages()) { print("not installed") }'], capture_output=True, text=True)

                            fh.write('#!/bin/bash\n')
                            fh.write(f'cd "/home/ayden/Documents/Silverman Lab/code/fenrir_paper_code/" \n\n')
                            if 'not installed' in result.stdout:
                                fh.write('Rscript -e "devtools::build(\'/home/ayden/Documents/Silverman Lab/code/fenrir/\',path=\'/home/ayden/Documents/Silverman Lab/code/fenrir/build/\')"\n')
                                fh.write('Rscript -e "install.packages(\'/home/ayden/Documents/Silverman Lab/code/fenrir/build/fenrir_1.0.tar.gz\', repos = NULL, type = \'source\')"\n')
                            
                            if infer_W == "0":
                                fh.write(f'Rscript code/cpp/optimizer.R {D_vals[D]} {Q_vals[Q]} "data/simulated/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}/" "results/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}_{m_idx}_{infer_W}_{run_mcmc}/" "1" {cores}\n')
                    
                            else:
                                fh.write(f'Rscript code/cpp/gibbs.R {D_vals[D]} {Q_vals[Q]} "data/simulated/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}/" "results/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}_{m_idx}_{infer_W}_{run_mcmc}/" {cores}\n')

                    call_str = f"sh {filename}"
                    print(f"Calling: {call_str}")
                    sys_response = subprocess.check_output(call_str, shell=True)
                    time.sleep(1)
                    sys_response = sys_response.decode().strip()
                    with open(f'job_listing_{vary}.txt', 'a') as job_listing:
                        job_listing.write(f"{sys_response[20:]}\tmodel={m_idx.upper()}\tN={N_vals[N]}\tD={D_vals[D]}\tQ={Q_vals[Q]}\tR={rseed[seed]}\n")


