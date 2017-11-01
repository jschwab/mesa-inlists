#!/bin/bash

cd ${SLURM_SUBMIT_DIR}
chmesadir r7624
./clean
./mk

# first make a pre-ms model
# ln -sf inlists/inlist_pms inlist_project
# ./rn

# then relax the composition
for f in inlists/inlist_composition_log*
do
   ln -sf inlists/inlist_composition inlist_project
   ln -sf $f inlist_composition_id
   ./rn
done

# and now go to the final density
ln -sf inlists/inlist_finalize_1em6 inlist_project
for f in models/composition_log*.mod
do
    FINAL=$(echo "$f" | sed 's/composition/final/')
    cp $f initial.mod
    ./rn
    cp final.mod $FINAL
done


