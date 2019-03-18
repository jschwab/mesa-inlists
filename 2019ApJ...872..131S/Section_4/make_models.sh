#!/bin/bash

export MESA_INLIST=inlist.simple

# first make a pre-ms model
ln -sf inlists/inlist_pms_tiny inlist_project
./rn

# then relax the composition
for f in inlists/inlist_composition_*
do
   ln -sf inlists/inlist_composition inlist_project
   ln -sf $f inlist_composition_id
   ./rn
done

# and now go to the final density
ln -sf inlists/inlist_finalize_1em6 inlist_project
for f in models/composition_*.mod
do
    FINAL=$(echo "$f" | sed 's/composition/final/')
    cp $f initial.mod
    ./rn
    cp final.mod $FINAL
done
