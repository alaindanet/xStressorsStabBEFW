#!/bin/bash


# Job name
#$ -N mod_comp_fw

# Error message
#$ -e mod_comp_fw.error
#$ -o mod_comp_fw.out

# Amount of RAM requested per node
#$ -l rmem=64G

# Replace by the path to the folder where your script lives if necessary
DIR_ENV=/home/${USER}/xStressorsStabBEFW
DIR_SCRIPT=script

# Load modules
module load apps/julia

cd ${DIR_ENV} && julia --project=${DIR_ENV} ${DIR_ENV}/${DIR_SCRIPT}/module.jl ${DIR_ENV}
