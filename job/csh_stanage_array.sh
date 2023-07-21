#!/bin/bash
#
#SBATCH --comment=CSZhK
#SBATCH --job-name=CSZhK
# Set number of iteration
#SBATCH --array=1-54%10
# Amount of RAM requested per job
#SBATCH --mem=32G
# Nb of threads requested per job (smp = shared memory)
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/users/${USER}/xStressorsStabBEFW
DIR_SCRIPT=scripts

# Load modules
#module load apps/julia/1.8.5/binary

# Put 1.8.5
JULIA=/users/bi1ahd/julia-1.9.2/bin/julia

STEPSIZE=4000

START=$((1 + ${STEPSIZE} * (${SLURM_ARRAY_TASK_ID} - 1)))
END=$((${STEPSIZE} * ${SLURM_ARRAY_TASK_ID}))

echo "Starting task from ${START} to ${END}"

cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/ct_S_h_stoch_array.jl ${START} ${END}
