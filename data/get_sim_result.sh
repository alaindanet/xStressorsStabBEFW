#!/bin/bash
#
SSH_ADDR=bi1ahd@stanage
FOLDER=/mnt/parscratch/users/bi1ahd/pre_process_xStressorsStabBEFW_targets/objects/

rsync -avP -e ssh --files-from=rsync-src-files ${SSH_ADDR}:${FOLDER} . #--dry-run
#rsync -avP -e ssh --files-from=rsync-src-files bi1ahd@stanage:/mnt/parscratch/users/bi1ahd/pre_process_xStressorsStabBEFW_targets/objects/ . --dry-run
