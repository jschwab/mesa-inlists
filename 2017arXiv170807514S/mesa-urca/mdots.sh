#!/bin/bash

export MESASDK_ROOT=/global/home/groups/ac_astro/jwschwab/mesasdk
source ${MESASDK_ROOT}/bin/mesasdk_init.sh
export MESA_DIR=/global/scratch/jwschwab/mesa-r9793

if [ -n "${SLURM_SUBMIT_DIR}" ]; then
   cd ${SLURM_SUBMIT_DIR}
fi

# rebuild
./clean
./mk

# use logft11 rxns
ln -sf transitions/weak-logft-11.transitions weak.transitions

# turn MLT off
ln -sf inlists/inlist_mlt_none inlist_variable

do_one_mdot() {

    # extract name
    NAME=mdot_${1}

    ln -sf inlists/inlist_${NAME} inlist_io    

    # remove old png files
    rm -rf png
    
    # run MESA
    ./rn

    # make movie
    images_to_movie.sh 'png/grid1*.png' ${NAME}-2017-06-23.mp4

}

do_one_mdot 3em6
do_one_mdot 1em6
do_one_mdot 3em7
do_one_mdot 1em7
do_one_mdot 3em8
do_one_mdot 1em8
do_one_mdot 3em9
do_one_mdot 1em9
