#!/bin/sh

for FILE in gasdens*.cfg ; do
  F=`echo $FILE | awk '{ print gensub("[a-z]*([0-9]*).cfg", "\\\\1", 1, $0); }'`
  echo $F
  ln -sf $FILE gasdens.cfg
  ln -sf gastemper$F.cfg gastemper.cfg
  ./viscosity.plt
  mv viscosity.png viscosity$F.png
done

qiv viscosity*.png

