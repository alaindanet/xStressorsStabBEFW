#!/bin/bash
#
#$ -N cs_big
#$ -t 1-277200:2000
#$ -o /home/$USER/log_test/
#$ -e /home/$USER/log_test/

# Load modules
module load apps/julia/1.8.5/binary


START=$SGE_TASK_ID
LAST=$(($SGE_TASK_ID + $SGE_TASK_STEPSIZE - 1))

echo "Starting task from ${START} to ${LAST}"
julia -e 'println("From $(ARGS[1]) to $(ARGS[2])")' ${START} ${LAST}

