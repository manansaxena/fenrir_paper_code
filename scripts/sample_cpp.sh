#!/bin/bash
#SBATCH -J mallard_cpp
#SBATCH --output=<random_path>/mallard_stan.out
#SBATCH --mem=256GB
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
module load anaconda
conda activate stanpip
module load r
cd "random_path_where_code_is_stored/fenrir_paper_code/" 

Rscript random_path_where_code_is_stored/fenrir_paper_code/code/cpp/mallard.R
