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

    ID=$(ruby make_inlist_io.rb ${NAME} ${OPTS})

    # remove old png files
    rm -rf png
    
    # run MESA
    ./rn

    # make movie
    images_to_movie.sh 'png/grid1*.png' movies/${ID}-2018-01-03.mp4

}

# at first, we don't need anything extra
ln -sf inlist_empty inlist_other

# stock case is mass loss and rotation, but not diffusion
ln -sf inlist_with_mass_loss inlist_mass_loss
ln -sf inlist_with_rotation inlist_rotation
ln -sf inlist_empty inlist_diffusion

do_one ZP6_XH_0.0_qH_0.75

do_one ZP6_XH_0.005_qH_0.75
do_one ZP6_XH_0.02_qH_0.75

do_one ZP6_XH_0.01_qH_0.6
do_one ZP6_XH_0.01_qH_0.7
do_one ZP6_XH_0.01_qH_0.75
do_one ZP6_XH_0.01_qH_0.8
do_one ZP6_XH_0.01_qH_0.9

do_one ZP7_XH_0.01_qH_0.75

# do two different mass loss cases
ln -sf inlist_with_alt_mass_loss inlist_mass_loss
ln -sf inlist_with_rotation inlist_rotation
ln -sf inlist_empty inlist_diffusion

do_one ZP6_XH_0.01_qH_0.75 alt
do_one ZP7_XH_0.01_qH_0.75 alt

ln -sf inlist_with_alt2_mass_loss inlist_mass_loss

do_one ZP6_XH_0.01_qH_0.75 alt2
do_one ZP7_XH_0.01_qH_0.75 alt2

# add diffusion
ln -sf inlist_with_mass_loss inlist_mass_loss
ln -sf inlist_with_rotation inlist_rotation
ln -sf inlist_with_diffusion inlist_diffusion

do_one ZP6_XH_0.01_qH_0.75 diff
do_one ZP7_XH_0.01_qH_0.75 diff

# without rotational mixing
ln -sf inlist_with_mass_loss inlist_mass_loss
ln -sf inlist_no_rotmix inlist_rotation
ln -sf inlist_empty inlist_diffusion

do_one ZP6_XH_0.01_qH_0.75 norotmix
do_one ZP7_XH_0.01_qH_0.75 norotmix

# without rotational mixing, with diffusion
ln -sf inlist_with_mass_loss inlist_mass_loss
ln -sf inlist_no_rotmix inlist_rotation
ln -sf inlist_with_diffusion inlist_diffusion

do_one ZP6_XH_0.01_qH_0.75 norotmix-diff
do_one ZP7_XH_0.01_qH_0.75 norotmix-diff

# without rotation or mass loss
ln -sf inlist_empty inlist_mass_loss
ln -sf inlist_no_rotation inlist_rotation
ln -sf inlist_empty inlist_diffusion

do_one ZP6_XH_0.01_qH_0.75 norot
do_one ZP7_XH_0.01_qH_0.75 norot

# without rotation or mass loss, with diffusion
ln -sf inlist_empty inlist_mass_loss
ln -sf inlist_no_rotation inlist_rotation
ln -sf inlist_with_diffusion inlist_diffusion

do_one ZP6_XH_0.01_qH_0.75 norot-diff


# stock case is mass loss and rotation, but not diffusion
ln -sf inlist_with_mass_loss inlist_mass_loss
ln -sf inlist_with_rotation inlist_rotation
ln -sf inlist_empty inlist_diffusion
ln -sf inlist_big_net inlist_other

do_one ZP6_XH_0.01_qH_0.75 big_net
