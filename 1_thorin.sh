#!/bin/bash

#cd src_main; ./make.sh; cd ..

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:src_reb

./thorin -vm in/in.par

#./thorin -vm -d in/in.par
#./thorin -vm -d -s 737 in/in_120.par
#./thorin -vm -d in/in_120.par

