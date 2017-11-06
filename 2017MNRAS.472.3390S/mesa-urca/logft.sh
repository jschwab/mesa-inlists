#!/bin/bash

source /global/home/users/jwschwab/mesa-init/mesasdk_init.sh
export MESA_DIR=/global/scratch/jwschwab/mesa-r9793

if [ -n "${SLURM_SUBMIT_DIR}" ]; then
   cd ${SLURM_SUBMIT_DIR}
fi
# rebuild
./clean
./mk

# do all this for SQB15+ model
ln -sf inlists/inlist_logft inlist_io
ln -sf inlists/inlist_mlt_none inlist_variable

# make a place to store the data
mkdir LOGS-logft

do_one_logft() {

    # extract name
    NAME=logft-${1}

    # link and run
    ln -sf transitions/weak-${NAME}.transitions weak.transitions
    ./rn

    # move data
    mv LOGS/history.data LOGS-logft/history-${NAME}.data
    mv LOGS/final_profile.data LOGS-logft/profile-${NAME}.data
    mv LOGS/final_model.mod LOGS-logft/model-${NAME}.mod

}

do_one_none() {

    # link and run
    ln -sf transitions/weak.transitions weak.transitions
    ./rn

    # move data
    NAME="no-forbidden"
    mv LOGS/history.data LOGS-logft/history-${NAME}.data
    mv LOGS/final_profile.data LOGS-logft/profile-${NAME}.data
    mv LOGS/final_model.mod LOGS-logft/model-${NAME}.mod

}


do_one_none
do_one_logft 11
do_one_logft 13
do_one_logft 15
