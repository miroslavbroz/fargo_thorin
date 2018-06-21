#!/bin/bash

cd src_main; ./make.sh; cd ..

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.

./thorin -vm -d in/in.par

