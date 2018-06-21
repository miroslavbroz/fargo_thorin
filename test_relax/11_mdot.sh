#!/bin/sh

for FILE in gasdens*.cfg ; do
  F=`echo $FILE | awk '{ print gensub("[a-z]*([0-9]*).cfg", "\\\\1", 1, $0); }'`
  echo $F
  ln -sf $FILE gasdens.cfg
  ln -sf gastemper$F.cfg gastemper.cfg
  ./mdot.plt
  mv mdot.png mdot$F.png
done

qiv mdot*.png

