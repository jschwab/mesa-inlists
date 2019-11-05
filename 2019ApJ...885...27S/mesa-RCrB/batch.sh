#!/bin/bash -x
#PBS -N RCrB
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -V
#PBS -m abe
#PBS -M jwschwab@ucsc.edu
#PBS -t 1-3

if [ -n "${SLURM_SUBMIT_DIR}" ]; then
   cd ${SLURM_SUBMIT_DIR}
   module load mesasdk/20180127
   export MESA_DIR=/pfs/jschwab/mesa-r11701
   export OMP_NUM_THREADS=16
else
    # assume you're already there
    SLURM_SUBMIT_DIR=$(pwd)
    SLURM_ARRAY_TASK_ID=1
fi

# rebuild MESA
# ./clean
#./mk

# copy AESOPUS opacity files
cp *.h5 ${MESA_DIR}/data/kap_data/

# make sure there's a place for the movies and output
mkdir -p movies
mkdir -p output

# define a function to run a named model
do_one() {

    # get the relevant line from the batch file
    read ID INLIST_MASS INLIST_Z INLIST_ABUNDANCES INLIST_OTHER <<< $(sed "${1}q;d" < $2)

    # use the main inlist in the submit directory
    export MESA_INLIST=${SLURM_SUBMIT_DIR}/inlist

    # make a temporary directory
    TMPDIR=$(mktemp -d)
    cd ${TMPDIR}

    # cache locally
    # mkdir -p caches
    # export MESA_CACHES_DIR=$(pwd)/caches

    # setup inlists via soft links
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_fixed .
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/${INLIST_MASS} inlist_mass
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/${INLIST_Z} inlist_metallicity
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_pgstar
    
    # softlink in some more stuff
    ln -sf ${SLURM_SUBMIT_DIR}/history_columns.list .
    ln -sf ${SLURM_SUBMIT_DIR}/profile_columns.list .
    ln -sf ${SLURM_SUBMIT_DIR}/zams .
    mkdir -p ${SLURM_SUBMIT_DIR}/output/${ID}
    ln -sf ${SLURM_SUBMIT_DIR}/output/${ID} LOGS

    # phase 1: HeMS to small He envelope
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/inlist_start inlist_variable
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/inlist_empty inlist_other
    
    # run MESA
    ${SLURM_SUBMIT_DIR}/star

    # phase 2: to the end!
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/inlist_rcb inlist_variable
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/${INLIST_ABUNDANCES} inlist_abundances
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/${INLIST_OTHER} inlist_other
    
    # restart from last photo
    PHOTO=$(ls -1 photos | tail -1)
    cp photos/${PHOTO} restart_photo
    ${SLURM_SUBMIT_DIR}/star

    # make movie
    DATE=$(date +%F)
    images_to_movie 'png/grid1*.png' ${SLURM_SUBMIT_DIR}/movies/${ID}-${DATE}.mp4

}

# do_one ${SLURM_ARRAY_TASK_ID} masses.batch
# do_one ${SLURM_ARRAY_TASK_ID} metallicity.batch
# do_one ${SLURM_ARRAY_TASK_ID} weiss.batch
# do_one ${SLURM_ARRAY_TASK_ID} bc.batch
# do_one ${SLURM_ARRAY_TASK_ID} A00.batch

