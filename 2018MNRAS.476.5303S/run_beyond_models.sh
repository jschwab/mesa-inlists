#!/bin/bash
#PBS -N sdb
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -V
#PBS -m abe
#PBS -M jwschwab@ucsc.edu

if [ -n "${PBS_O_WORKDIR}" ]; then
    cd ${PBS_O_WORKDIR}
    module load mesasdk/20160129
    export MESA_DIR=/pfs/jschwab/mesa-r10108
fi

# rebuild MESA
./clean
./mk

# make sure there's a place for the movies
mkdir -p movies

# define a function to run a named model
do_one() {

    # extract name
    NAME=${1}
    OPTS=${2}

    ID=$(ruby make_beyond_inlist_io.rb ${NAME} ${OPTS})

    # remove old png files
    rm -rf png
    
    # run MESA
    ./rn

    # make movie
    images_to_movie.sh 'png/grid1*.png' movies/${ID}-2018-01-03.mp4

}

# go to WD phase
ln -sf inlist_beyond inlist_other

# stock case is rotation, no mass loss, and no diffusion
ln -sf inlist_empty inlist_mass_loss
ln -sf inlist_with_rotation inlist_rotation
ln -sf inlist_empty inlist_diffusion

do_one ZP6_XH_0.01_qH_0.75-sdb WD
do_one ZP7_XH_0.01_qH_0.75-sdb WD