#!/bin/bash

source /global/home/users/jwschwab/mesa-init/mesasdk_init.sh
export MESA_DIR=/global/scratch/jwschwab/mesa-MR16

./clean
./mk

ln -sf inlists/inlist_simmer inlist_fixed

do_one() {

    # link relevant inlists
    ln -sf inlists/inlist_${1} inlist_variable

    # do run
    ./rn

    # save logs
    rm -rf LOGS-${2}
    mv LOGS LOGS-${2}
    
    # make movie
    images_to_movie.sh 'png/grid1*.png' ${2}-2017-06-15.mp4
    rm -rf png

}

ln -sf inlists/inlist_io_toy inlist_io
do_one MR16_rates MR16
do_one MR16_mixing MR16-adv
do_one MR16_mixing_convergence MR16-adv-hr

ln -sf inlists/inlist_io_toy-lowC inlist_io
do_one MR16_rates MR16-lowC
