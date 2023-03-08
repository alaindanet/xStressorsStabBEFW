#!/bin/bash
#
#$ -N test
#$ -t 1-1000:500
#$ -o /home/$USER/log_test/
#$ -e /home/$USER/log_test/

# Load modules
module load apps/julia/1.8.5/binary


START=$SGE_TASK_ID
LAST=$(($SGE_TASK_ID + $SGE_TASK_STEPSIZE - 1))

echo "Starting task from ${START} to ${LAST}"
julia -e 'println("From $(ARGS[1]) to $(ARGS[2])")' $SGE_TASK_ID $(($SGE_TASK_ID + $SGE_TASK_STEPSIZE - 1))

