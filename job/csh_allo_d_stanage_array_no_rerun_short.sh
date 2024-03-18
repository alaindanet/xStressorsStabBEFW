#!/bin/bash
#
#SBATCH --job-name=allo_d_no_rerun_short
# Set number of iteration
#SBATCH --array=1-34%20
# Amount of RAM requested per job
#SBATCH --mem=64G
# Nb of threads requested per job (smp = shared memory)
#SBATCH --cpus-per-task=20
#SBATCH --ntasks=1
# Request
#SBATCH --time=2-00:00:00

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/users/bi1ahd/xStressorsStabBEFW
DIR_SCRIPT=scripts

# Load modules
#module load apps/julia/1.8.5/binary

# Put 1.8.5
JULIA="/users/bi1ahd/.juliaup/bin/julialauncher +release"

STEPSIZE=2000

START=$((1 + ${STEPSIZE} * (${SLURM_ARRAY_TASK_ID} - 1)))
END=$((${STEPSIZE} * ${SLURM_ARRAY_TASK_ID}))

echo "Starting task from ${START} to ${END}"

cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/ct_S_h_d_allometric_array_no_rerun_short.jl ${START} ${END}
