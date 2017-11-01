#!/bin/bash

source /global/home/users/jwschwab/mesa-init/mesasdk_init.sh
export MESA_DIR=/global/scratch/jwschwab/mesa-r9793

if [ -n "${SLURM_SUBMIT_DIR}" ]; then
   cd ${SLURM_SUBMIT_DIR}
fi
# rebuild
./clean
./mk

# define a function to run a named model
do_one_named() {

    # extact name
    NAME=${1}

    ln -sf inlists/inlist_io_${NAME} inlist_io
    ln -sf inlists/inlist_${NAME} inlist_variable
    ./rn

}

# use logft11 rxns
ln -sf transitions/weak-logft-11.transitions weak.transitions

#do_one_named SQB15+
#do_one_named T13
#do_one_named F15

# use no forbidden rxns
ln -sf transitions/weak.transitions weak.transitions

do_one_named SQB15
