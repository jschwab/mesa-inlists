#!/bin/bash

#PBS -N composition
#PBS -l nodes=1:ppn=16
#PBS -l walltime=48:00:00
#PBS -V
#PBS -j oe
#PBS -m abe
#PBS -M jwschwab@ucsc.edu


# setup hyades specific stuff
if [ -n "${PBS_O_WORKDIR}" ]; then
    # load SDK
    module load mesasdk/20180127
    export MESA_DIR=/pfs/jschwab/mesa-r10108
    cd ${PBS_O_WORKDIR}
fi

# rebuild MESA
./clean
./mk

# define a function to run a named model
do_one_composition() {

    # extact name
    NAME=${1}

    ln -sf inlists/inlist_io_composition_${NAME} inlist_io
    ./rn

}

# use logft11 rxns
ln -sf transitions/weak-logft-11.transitions weak.transitions

# turn off MLT
ln -sf inlists/inlist_mlt_none inlist_variable


do_one_composition nourca_xc-0.000
do_one_composition nourca_xc-0.001
do_one_composition nourca_xc-0.003
do_one_composition nourca_xc-0.010
do_one_composition nourca_xc-0.020
do_one_composition nourca_xc-0.030
do_one_composition nourca_xc-0.040
do_one_composition nourca_xc-0.050
do_one_composition nourca_xc-0.100


do_one_composition withurca_xc-0.000
do_one_composition withurca_xc-0.001
do_one_composition withurca_xc-0.003
do_one_composition withurca_xc-0.010
do_one_composition withurca_xc-0.020
do_one_composition withurca_xc-0.030
do_one_composition withurca_xc-0.040
do_one_composition withurca_xc-0.050
do_one_composition withurca_xc-0.100
