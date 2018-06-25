#!/bin/bash

#cd src_main; ./make.sh; cd ..

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:src_reb

./thorin -vm in_relax/in.par


