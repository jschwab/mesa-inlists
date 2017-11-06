#!/bin/bash

source /global/home/users/jwschwab/mesa-init/mesasdk_init.sh
export MESA_DIR=/global/scratch/jwschwab/mesa-r9793

if [ -n "${SLURM_SUBMIT_DIR}" ]; then
   cd ${SLURM_SUBMIT_DIR}
fi

./clean
./mk

for f in inlists-scalings/inlist_*
do
   ln -sf $f inlist_variable
   ./rn
done
