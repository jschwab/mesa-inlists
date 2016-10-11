#!/bin/bash

# MESA version r6596

# recompile MESA
./clean
./mk

# evolve a bunch of models

ln -sf inlists/inlist_evolve inlist_fixed

for f in inlists/inlist_evolve_*
do
   ln -sf $f inlist_variable
   ./rn
done

