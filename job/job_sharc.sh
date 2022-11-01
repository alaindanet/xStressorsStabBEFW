#!/bin/bash                                                                                                               

# Job name 
#$ -N vasseur_fox_2007 

# Error message
#$ -e vasseur_fox_2007.error
#$ -o vasseur_fox_2007.out

# Amount of RAM requested per node 
#$ -l rmem=16G

# Replace by the path to the folder where your script lives if necessary  
DIR_SCRIPT=/home/${USER}/xStressorsStabBEFW

# Load modules 
module load apps/julia

julia --project=${DIR_SCRIPT} ${DIR_SCRIPT}/vasseur_fox_2007.jl
