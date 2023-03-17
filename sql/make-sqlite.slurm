#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE        #Do not propagate environment
#SBATCH --get-user-env=L     #Replicate login environment
  
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=sql-metaT     #Set the job name to "JobExample1"
#SBATCH --time=08:00:00            #Set the wall clock limit to 1hr and 30min
#SBATCH --ntasks=4                 #Request 1 task
#SBATCH --ntasks-per-node=64        #Request 1 task/core per node
#SBATCH --mem=248G                #Request 2560MB (2.5GB) per node,250G
#SBATCH --output=sql.%j    #Send stdout/err to "Example1Out.[jobID]"
#SBATCH --partition=cpu

export OMP_NUM_THREADS=16

# To use conda, include this line:
eval "$(conda shell.bash hook)"

echo "load modules"
module load GCCcore/12.2.0
module load SQLite/3.39.4
module load Python/3.10.8

module load GCC/11.2.0
module load OpenMPI/4.1.1
module load pandas/1.4.0-Python-3.9.6


echo "launch python"

python make_sql_db-toy.py

echo "completed adding tables to sqlite database"

