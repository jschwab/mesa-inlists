#!/bin/bash -x
#SBATCH --job-name=RCrB
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH --time=08:00:00
#SBATCH --array 1-51
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
    SLURM_ARRAY_TASK_ID=19
fi

# rebuild MESA
#./clean
#./mk

# make sure there's a place for the output and final models
mkdir -p models
mkdir -p output
mkdir -p final_models

# define a function to run a named model
do_one() {

    # get the relevant line from the batch file
    read ID <<< $(sed "${1}q;d" < $2)

    # use the main inlist in the submit directory
    export MESA_INLIST=${SLURM_SUBMIT_DIR}/inlist

    # make a temporary directory
    TMPDIR=$(mktemp -d)
    cd ${TMPDIR}

    # cache locally
    # mkdir -p caches
    # export MESA_CACHES_DIR=$(pwd)/caches

    # setup inlists via soft links
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_common .
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_pgstar .
    
    # softlink in some more stuff
    ln -sf ${SLURM_SUBMIT_DIR}/compositions .
    ln -sf ${SLURM_SUBMIT_DIR}/models .
    mkdir -p ${SLURM_SUBMIT_DIR}/output/${ID}
    ln -sf ${SLURM_SUBMIT_DIR}/output/${ID} LOGS
    ln -sf ${SLURM_SUBMIT_DIR}/final_models .

    # phase 1: make a pre-ms model
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_pms inlist_fixed
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/inlist_pms_${ID} inlist_variable
    
    # run MESA
    ${SLURM_SUBMIT_DIR}/star

    # phase 2: relax composition
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_composition inlist_fixed
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/inlist_composition_${ID} inlist_variable
    
    # run MESA
    ${SLURM_SUBMIT_DIR}/star

    # phase 3: relax entropy
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_entropy inlist_fixed
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/inlist_entropy_${ID} inlist_variable

    # run MESA
    ${SLURM_SUBMIT_DIR}/star

    # phase 4: finalize model
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_finalize inlist_fixed
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/inlist_finalize_${ID} inlist_variable

    # run MESA
    ${SLURM_SUBMIT_DIR}/star

    cd -
    rm -rf ${TMPDIR}
    
}

do_one ${SLURM_ARRAY_TASK_ID} models.batch
