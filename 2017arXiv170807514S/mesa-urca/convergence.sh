#!/bin/bash

source /global/home/users/jwschwab/mesa-init/mesasdk_init.sh
export MESA_DIR=/global/scratch/jwschwab/mesa-r9793

if [ -n "${SLURM_SUBMIT_DIR}" ]; then
   cd ${SLURM_SUBMIT_DIR}
fi

./clean
./mk
ln -sf inlists/inlist_empty inlist_io

# loop over inlists
for f in inlists/inlist_convergence_*
do
   ln -sf $f inlist_variable
   ./rn
done
