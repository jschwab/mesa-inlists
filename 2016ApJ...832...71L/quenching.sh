#!/bin/bash

./clean
./mk

do_one() {

    # link inlists
    ln -sf inlists/inlist_${1} inlist_variable
    
    # run MESA
    ./rn

    # save history
    cp LOGS/history.data ${1}-history.data
    
    # make movie
    images_to_movie.sh 'png/profile*.png' ${1}-profile.mp4
    images_to_movie.sh 'png/History*.png' ${1}-history.mp4
    rm png/*.png

}

do_one ncrit_0.001
do_one ncrit_0.003
do_one ncrit_0.010
do_one ncrit_0.030
do_one ncrit_0.100
do_one ncrit_0.300
do_one ncrit_0.500

do_one Le-0.1_Nc-0.100
do_one Le-0.1_Nc-0.150
do_one Le-0.1_Nc-0.200
do_one Le-0.1_Nc-0.250
do_one Le-0.3_Nc-0.100
do_one Le-0.3_Nc-0.150
do_one Le-0.3_Nc-0.200
do_one Le-0.3_Nc-0.250
do_one Le-0.5_Nc-0.100
do_one Le-0.5_Nc-0.150
do_one Le-0.5_Nc-0.200
do_one Le-0.5_Nc-0.250
do_one Le-1.0_Nc-0.100
do_one Le-1.0_Nc-0.150
do_one Le-1.0_Nc-0.200
do_one Le-1.0_Nc-0.250
do_one Le-3.0_Nc-0.100
do_one Le-3.0_Nc-0.150
do_one Le-3.0_Nc-0.200
do_one Le-3.0_Nc-0.250
do_one Le-5.0_Nc-0.100
do_one Le-5.0_Nc-0.150
do_one Le-5.0_Nc-0.200
do_one Le-5.0_Nc-0.250
do_one Le-10.0_Nc-0.100
do_one Le-10.0_Nc-0.150
do_one Le-10.0_Nc-0.200
do_one Le-10.0_Nc-0.250
