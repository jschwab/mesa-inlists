#!/bin/bash

# MESA version r6794

# recompile MESA
./clean
./mk

# make a bunch of pre-ms models

ln -sf inlists/inlist_pms inlist_fixed

for f in inlists/inlist_pms_1p*
do
   ln -sf $f inlist_variable
   ./rn
done

