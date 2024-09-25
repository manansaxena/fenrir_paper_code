#!/bin/bash
#SBATCH -J N300_D3_Q1_R1_mycpp
#SBATCH --output=<random_path>/N300_D3_Q1_R1_mycpp.out
#SBATCH --mem=256GB
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
module load anaconda
conda activate stanpip
module load r
cd "random_path_where_code_is_stored/fenrir_paper_code/" 

Rscript random_path_where_code_is_stored/fenrir_paper_code/code/cpp/optimizer.R 3 1 "random_path_where_code_is_stored/fenrir_paper_code/data/simulated/N300_D3_Q1_R1_W0.45_pm5/" "random_path_where_results_are_stored/N300_D3_Q1_R1_W0.45_pm5_mycpp_0_1/" "0" 16
