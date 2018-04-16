#!/bin/bash

# loop and check each model
for f in zeus_models/*_DT.dat
do
    # extract ID from file name
    ID=$(basename $f _DT.dat)
    if [ ! -f models/${ID}.mod ]; then
        echo ${ID} failed
    fi
done
