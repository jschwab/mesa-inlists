#!/bin/bash

export MESA_INLIST=inlist_make_model

do_one() {

    # get the run id
    ID=${1}

    # copy the 1D profiles
    cp zeus_models/${ID}_j.dat zeus_j.dat
    cp zeus_models/${ID}_X.dat zeus_X.dat
    cp zeus_models/${ID}_DT.dat zeus_DT.dat
    cp zeus_models/${ID}_s.dat zeus_s.dat

    # have MESA do the relaxation
    ./rn
    
    # save output
    mv final.mod models/${ID}.mod
    mv final.profile models/${ID}.profile

    # clean up
    rm zeus_j.dat zeus_X.dat zeus_DT.dat    
}


ln -sf initial_mass_ZP6 initial_mass

# loop and make each model
for f in zeus_models/ZP6_XH*_DT.dat
do
    # extract ID from file name
    ID=$(basename $f _DT.dat)
    do_one ${ID}
done


ln -sf initial_mass_ZP7 initial_mass

# loop and make each model
for f in zeus_models/ZP7_XH*_DT.dat
do
    # extract ID from file name
    ID=$(basename $f _DT.dat)
    do_one ${ID}
done

