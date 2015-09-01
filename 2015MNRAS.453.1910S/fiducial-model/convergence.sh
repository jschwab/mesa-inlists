#!/bin/bash

do_one () {
    ln -sf inlists/inlist_$1 inlist_variable
    ./rn
}

do_one both
do_one fiducial
do_one spatial
do_one temporal
