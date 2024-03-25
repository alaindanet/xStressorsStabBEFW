#!/bin/bash
#
#SBATCH --job-name=allo_d
# Set number of iteration
#SBATCH --array=1-29%20
# Amount of RAM requested per job
#SBATCH --mem=64G
# Nb of threads requested per job (smp = shared memory)
#SBATCH --cpus-per-task=20
#SBATCH --ntasks=1
# Request
#SBATCH --time=2-00:00:00

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/users/${USER}/xStressorsStabBEFW
DIR_SCRIPT=scripts

# Load modules
#module load apps/julia/1.8.5/binary

# Put 1.8.5
JULIA="/users/bi1ahd/.juliaup/bin/julialauncher +release"

STEPSIZE=2000

START=$((1 + ${STEPSIZE} * (${SLURM_ARRAY_TASK_ID} - 1)))
END=$((${STEPSIZE} * ${SLURM_ARRAY_TASK_ID}))

echo "Starting task from ${START} to ${END}"

cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/simulation_args.jl \
    --first_sim=${START} --last_sim=${END}\
    --param_file="scripts/param_comb_zc.arrow"\
    --save_dir="/mnt/parscratch/users/bi1ahd/sim/sim_csc_allo_d1/"\
    --tmax=2000\
    --K_corrected=true --K=10.0\
    --h=2.0\
    --rebuild_after_disconnected=false\
    --re_run=true\
    --d=nothing\
    --d_allometric_set="(ap = 0.4, ai = 0.4, ae = 0.4)"

cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/simulation_args.jl \
    --first_sim=${START} --last_sim=${END}\
    --param_file="scripts/param_comb_ct_S_h_d4.arrow"\
    --save_dir="/mnt/parscratch/users/bi1ahd/sim/sim_csc_allo_d1_K_no_corrected/"\
    --tmax=2000\
    --K_corrected=false --K=1.0\
    --h=2.0\
    --rebuild_after_disconnected=false\
    --re_run=true\
    --d=nothing\
    --d_allometric_set="(ap = 0.4, ai = 0.4, ae = 0.4)"

cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/simulation_args.jl \
    --first_sim=${START} --last_sim=${END}\
    --param_file="scripts/param_comb_zc.arrow"\
    --save_dir="/mnt/parscratch/users/bi1ahd/sim/sim_csc_allo_d1_non_allo/"\
    --tmax=2000\
    --K_corrected=true --K=10.0\
    --h=2.0\
    --rebuild_after_disconnected=false\
    --re_run=true\
    --d=.1\
    --d_allometric_set=nothing
