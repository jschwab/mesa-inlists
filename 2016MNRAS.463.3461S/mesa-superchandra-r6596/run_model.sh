#!/bin/bash

# MESA version r6596

./clean
./mk

do_one() {

    # make temporary working directory
    WORKDIR=$(mktemp -d --tmpdir=.)
    cd $WORKDIR

    # copy inlists
    cp ../inlist inlist
    cp ../inlists/inlist_${1} inlist_fixed
    cp ../inlists/inlist_${1}-${2} inlist_variable
    cp ../inlist_pgstar .
    
    # copy column files
    cp ../history_columns.list .
    cp ../profile_columns.list .
    
    # copy executable files
    cp ../re .
    cp ../rn .
    cp ../star .

    # copy other files
    ln -sf ../models .
    cp ../CO-table-with-metals.dat .

    # link output dir
    ln -sf ../mesa_output .

    # rn
    ./rn

    cd ..

}

do_one cflame M15
do_one kh M15
do_one neflame M15
do_one khNe M15
