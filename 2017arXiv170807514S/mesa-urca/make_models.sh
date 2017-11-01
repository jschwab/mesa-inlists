#!/bin/bash

export MESA_INLIST=inlist.simple

# first make a pre-ms model
ln -sf inlists/inlist_pms inlist_project
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


# do the same procedure for a low mass model
ln -sf inlists/inlist_pms_small inlist_project
./rn

ln -sf inlists/inlist_composition inlist_project
ln -sf inlists/inlist_composition_SQB15+ inlist_composition_id
./rn

ln -sf inlists/inlist_finalize_3em8 inlist_project
cp models/composition_SQB15+.mod initial.mod
./rn
cp final.mod models/final_SQB15+_small_3em8.mod
