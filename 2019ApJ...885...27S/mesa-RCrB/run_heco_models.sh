#!/bin/bash -x
#SBATCH --job-name=RCrB
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH --time=08:00:00
#SBATCH --array 1-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jwschwab@ucsc.edu
#SBATCH -A csc116

if [ -n "${SLURM_SUBMIT_DIR}" ]; then
   cd ${SLURM_SUBMIT_DIR}
   module load mesasdk/20190503
   export MESA_DIR=${PROJECT}/mesa-r11701
   export MESA_CACHES_DIR=/oasis/scratch/comet/jschwab/temp_project/mesa-r11701/cache
   export OMP_NUM_THREADS=4
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
    read ID INLIST_LOAD INLIST_Z INLIST_OTHER <<< $(sed "${1}q;d" < $2)

    # use the main inlist in the submit directory
    export MESA_INLIST=${SLURM_SUBMIT_DIR}/inlist_heco

    # make a temporary directory
    TMPDIR=$(mktemp -d)
    cd ${TMPDIR}

    # cache locally
    # mkdir -p caches
    # export MESA_CACHES_DIR=$(pwd)/caches

    # setup inlists via soft links
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_fixed .
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/${INLIST_Z} inlist_metallicity
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_pgstar .
    
    # softlink in some more stuff
    ln -sf ${SLURM_SUBMIT_DIR}/history_columns.list .
    ln -sf ${SLURM_SUBMIT_DIR}/profile_columns.list .
    ln -sf ${SLURM_SUBMIT_DIR}/models .
    ln -sf ${SLURM_SUBMIT_DIR}/heco_models .
    mkdir -p ${SLURM_SUBMIT_DIR}/output/${ID}
    ln -sf ${SLURM_SUBMIT_DIR}/output/${ID} LOGS

    # load saved model
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/inlist_evolve_heco inlist_variable
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/${INLIST_LOAD} inlist_load
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/${INLIST_OTHER} inlist_other
    
    # restart from last photo
    ${SLURM_SUBMIT_DIR}/star

    # make movie
    DATE=$(date +%F)
    images_to_movie 'png/grid1*.png' ${SLURM_SUBMIT_DIR}/movies/${ID}-${DATE}.mp4

    cd -
    rm -rf ${TMPDIR}
    
}

#do_one ${SLURM_ARRAY_TASK_ID} ZP4.batch
do_one ${SLURM_ARRAY_TASK_ID} heco.batch


