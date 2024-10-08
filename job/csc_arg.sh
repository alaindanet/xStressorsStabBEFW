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

# Default simulation
cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/simulation_args.jl \
    --first_sim=${START} --last_sim=${END}\
    --param_file="scripts/param_comb_zc.arrow"\
    --save_dir="/mnt/parscratch/users/bi1ahd/sim/sim_csc_allo_d2/"\
    --tmax=2000\
    --K_corrected=true --K=10.0\
    --h=2.0\
    --rebuild_after_disconnected=false\
    --re_run=true\
    --d=nothing

# Rebuild food-web after removing disconnected species
cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/simulation_args.jl \
    --first_sim=${START} --last_sim=${END}\
    --param_file="scripts/param_comb_zc.arrow"\
    --save_dir="/mnt/parscratch/users/bi1ahd/sim/sim_csc_allo_d2_rebuild/"\
    --tmax=2000\
    --K_corrected=true --K=10.0\
    --h=2.0\
    --rebuild_after_disconnected=true\
    --re_run=true\
    --d=nothing

# Do not standardise carrying capacity according the number of primary producers
cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/simulation_args.jl \
    --first_sim=${START} --last_sim=${END}\
    --param_file="scripts/param_comb_zc.arrow"\
    --save_dir="/mnt/parscratch/users/bi1ahd/sim/sim_csc_allo_d2_K_no_corrected/"\
    --tmax=2000\
    --K_corrected=false --K=10.0\
    --h=2.0\
    --rebuild_after_disconnected=false\
    --re_run=true\
    --d=nothing

# Death rates are fixed, i.e. are not allometric
cd ${DIR_ENV} && ${JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/simulation_args.jl \
    --first_sim=${START} --last_sim=${END}\
    --param_file="scripts/param_comb_zc.arrow"\
    --save_dir="/mnt/parscratch/users/bi1ahd/sim/sim_csc_allo_d2_non_allo/"\
    --tmax=2000\
    --K_corrected=true --K=10.0\
    --h=2.0\
    --rebuild_after_disconnected=false\
    --re_run=true\
    --d=0.1
