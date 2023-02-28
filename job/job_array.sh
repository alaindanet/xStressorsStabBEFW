#!/bin/bash
#
#$ -N Test array
#$ -t 1-1000:100
#$ -o /home/$USER/logs/
#$ -e /home/$USER/logs/

# Load modules
module load apps/julia/1.8.5/binary


echo "Starting task $SGE_TASK_ID"
for subtask_id in $(seq $SGE_TASK_ID $(( $SGE_TASK_ID + $SGE_TASK_STEPSIZE - 1 )) ); do
    echo "Subtask $subtask_id of task $SGE_TASK_ID"
    julia -e 'println($ARGS)' $subtask_id
done
