#!/bin/bash
#
#$ -N CSZhK
#
# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/home/${USER}/xStressorsStabBEFW
DIR_SCRIPT=scripts

# Set number of iteration
#IDNB=($(wc -l ${DIR_ENV}/${DIR_SCRIPT}/param_comb_ct_S_h.csv))
#IDNB=
#$ -t 1-288000:4000
#$ -o /home/$USER/logs/
#$ -e /home/$USER/logs/
#
# Amount of RAM requested per job
#$ -l rmem=64G
# Nb of threads requested per job (smp = shared memory)
#$ -pe smp 15


# Load modules
module load apps/julia/1.8.5/binary

# Put 1.8.5
#MY_JULIA=/home/bi1ahd/julia-1.8.5/bin/julia
MY_JULIA=julia


START=$SGE_TASK_ID
END=$(($SGE_TASK_ID + $SGE_TASK_STEPSIZE - 1))

echo "Starting task from ${START} to ${END}"
cd ${DIR_ENV} && ${MY_JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/ct_S_h_stoch_array.jl ${START} ${END}
