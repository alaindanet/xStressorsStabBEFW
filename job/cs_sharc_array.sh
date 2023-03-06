#!/bin/bash
#
#$ -N myArrayTest
#$ -t 1-18000:2000
#$ -o /home/$USER/logs/
#$ -e /home/$USER/logs/
#
# Amount of RAM requested per node
#$ -l rmem=64G

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/home/${USER}/xStressorsStabBEFW
DIR_SCRIPT=scripts

# Load modules
module load apps/julia/1.8.5/binary

# Put 1.8.5
#MY_JULIA=/home/bi1ahd/julia-1.8.5/bin/julia
MY_JULIA=julia




echo "Starting task from $SGE_TASK_ID to $($SGE_TASK_ID + $SGE_TASK_STEPSIZE - 1)"
cd ${DIR_ENV} && ${MY_JULIA} --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/connectance_richness_stoch_array.jl $SGE_TASK_ID $($SGE_TASK_ID + $SGE_TASK_STEPSIZE - 1)
