#!/bin/bash

# Job name
#$ -N mccann_brose

# Error message
#$ -e mccann_brose.error
#$ -o mccann_brose.out

# Amount of RAM requested per node
#$ -l rmem=16G

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/home/${USER}/xStressorsStabBEFW
DIR_SCRIPT=use_case

# Load modules
module load apps/julia

julia --project=${DIR_SCRIPT} ${DIR_ENV}/${DIR_SCRIPT}/vasseur_fox_brose.jl #vasseur_fox_2007.jl
