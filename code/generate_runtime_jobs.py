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
batch_jobs_dir = 'batch_jobs'
code_dir = "random_path_where_code_is_stored"
data_dir = "random_path_where_data_is_stored"
results_dir = "random_path_where_results_are_stored"
scripts_dir = "random_path_where_scripts_are_stored"

if vary == 'N':
    N_vals = ['200','500', '1000', '2000', '4000']
    D_vals = ['30']
    Q_vals = ['1']
    
elif vary == 'D':
    N_vals = ['600']
    D_vals = ['25', '50', '75', '100']
    Q_vals = ['1']

N_total_list, M0_list, C0_list, num_timeseries = generate_params(N_vals)
W_val = "0.45"
percent_of_missing = "5"
infer_W = "0"
run_mcmc = "0"

methods = ['sc','mycpp']
rseed = ['1','2','3','4','5','6','7','8','9','10']

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
                simulation_ret = subprocess.call(["Rscript", f"{code_dir}/fenrir_paper_code/code/generate_simulated_data.R", D_vals[D], Q_vals[Q], rseed[seed], W_val, percent_of_missing, N_total_str, M0_str, C0_str, f"{data_dir}/simulated"])
                if simulation_ret != 0:
                    print(f"Problem generating simulated data for {N_vals[N]}, {D_vals[D]}, {Q_vals[Q]} case (varying {vary})!")
                    exit(1)

                for m_idx in methods:

                    filename = f'{scripts_dir}/scripts/{vary}-varying_N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}_{m_idx}_{infer_W}_{run_mcmc}.sh'
                    
                    with open(filename, 'w') as fh:
                        if not os.path.exists(f'{results_dir}/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}_{m_idx}_{infer_W}_{run_mcmc}'):
                                os.mkdir(f'{results_dir}/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}_{m_idx}_{infer_W}_{run_mcmc}')

                        fh.write('#!/bin/bash\n')
                        fh.write(f'#SBATCH -J N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_{m_idx}\n')
                        fh.write(f'#SBATCH --output={batch_jobs_dir}/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_{m_idx}.out\n')
                        fh.write(f'#SBATCH --mem=256GB\n')
                        fh.write(f'#SBATCH --nodes=1\n')
                        fh.write(f'#SBATCH --ntasks={cores}\n')
                        fh.write(f'#SBATCH --time=48:00:00\n')
                        fh.write(f'module load anaconda\n')
                        fh.write(f'conda activate stanpip\n')
                        fh.write(f'module load r\n')

                        fh.write(f'cd {code_dir} \n\n')
                        
                        if m_idx == "sc":
                            fh.write(f'random_path_where_python_is_stored/python {code_dir}/fenrir_paper_code/code/stan/run_stan.py {D_vals[D]} {Q_vals[Q]} "{data_dir}/simulated/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}/" "{results_dir}/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}_{m_idx}_{infer_W}_{run_mcmc}/" "{code_dir}/fenrir_paper_code/code/stan/" {infer_W} {run_mcmc}')

                        elif m_idx == "mycpp": 

                            result = subprocess.run(['Rscript', '-e', 'if (!"fenrir" %in% installed.packages()) { print("not installed") }'], capture_output=True, text=True)
                            if 'not installed' in result.stdout:
                                fh.write('Rscript -e "devtools::build(\'random_path_where_fenrir_repo_is_stored/fenrir/\',path=\'random_path_where_fenrir_repo_is_stored/fenrir/build/\')"\n')
                                fh.write('Rscript -e "install.packages(\'random_path_where_fenrir_repo_is_stored/fenrir/build/fenrir_1.0.tar.gz\',  lib=\'random_path_where_R_libs_are_stored/R/x86_64-pc-linux-gnu-library/4.2\',  repos = NULL, type = \'source\')"\n')
                            
                            if infer_W == "0":
                                fh.write(f'Rscript {code_dir}/fenrir_paper_code/code/cpp/optimizer.R {D_vals[D]} {Q_vals[Q]} "{data_dir}/simulated/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}/" "{results_dir}/N{N_vals[N]}_D{D_vals[D]}_Q{Q_vals[Q]}_R{rseed[seed]}_W{W_val}_pm{percent_of_missing}_{m_idx}_{infer_W}_{run_mcmc}/" {cores}\n')

                    call_str = f"sbatch {filename}"
                    print(f"Calling: {call_str}")
                    sys_response = subprocess.check_output(call_str, shell=True)
                    time.sleep(1)
                    sys_response = sys_response.decode().strip()
                    with open(f'job_listing_temp_{vary}.txt', 'a') as job_listing:
                        job_listing.write(f"{sys_response[20:]}\tmodel={m_idx.upper()}\tN={N_vals[N]}\tD={D_vals[D]}\tQ={Q_vals[Q]}\tR={rseed[seed]}\n")


