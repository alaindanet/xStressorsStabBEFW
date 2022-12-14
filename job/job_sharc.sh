#!/bin/bash


# Job name
#$ -N cs_stoch 

# Error message
#$ -e cs_stoch.error
#$ -o cs_stoch.out

# Amount of RAM requested per node
#$ -l rmem=32G

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/home/${USER}/xStressorsStabBEFW
DIR_SCRIPT=use_case

# Load modules
module load apps/julia

julia --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/connectance_richness_stoch.jl ${DIR_ENV} #vasseur_fox_2007.jl
