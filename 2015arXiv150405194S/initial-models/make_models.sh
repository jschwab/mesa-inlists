#!/bin/bash

# first make a pre-ms model
ln -sf inlists/inlist_pms inlist_project
./rn

# then relax the composition
for f in inlists/inlist_composition_??????
do
   ln -sf $f inlist_project
   ./rn
done

# and now go to final density, but with a variety of mdots, so that
# you get a variety of central temperatures

for f in inlists/inlist_finalize_1em?
do
    ln -sf $f inlist_project
    MDOT=$(echo -n "$f" | tail -c 4)
    mkdir -p models/mdot_$MDOT

    for g in models/composition_??????.mod
    do
        FINAL=$(echo "$g" | sed 's/composition/final/')
        cp $g initial.mod
        ./rn
        cp final.mod $FINAL
        cp $FINAL models/mdot_$MDOT/
    done

done

