#!/bin/bash
#
#$ -N update 
#$ -o /home/$USER/log_test/
#$ -e /home/$USER/log_test/

# Amount of RAM requested per node
#$ -l rmem=64G

# Load modules
module load apps/julia/1.8.5/binary



cd /home/$USER/BEFWM2_fork/ && julia --procs=1 --project=. -e 'import Pkg; Pkg.update()'
cd /home/$USER/xStressorsStabBEFW/ && julia --procs=1 --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.update()'
cd /home/$USER/xStressorsStabBEFW/ && julia --procs=1 --project=. -e 'using EcologicalNetworksDynamics; FoodWeb(nichemodel, 10, C = .1)'

