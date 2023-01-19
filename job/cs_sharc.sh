#!/bin/bash


# Job name
#$ -N cs

# Error message
#$ -e cs.error
#$ -o cs.out

# Amount of RAM requested per node
#$ -l rmem=64G

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/home/${USER}/xStressorsStabBEFW
DIR_SCRIPT=scripts

# Load modules
module load apps/julia

cd ${DIR_ENV} && julia --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/connectance_richness_stoch.jl ${DIR_ENV}
