#!/bin/sh

for FILE in gastemper*.cfg ; do
  F=`echo $FILE | awk '{ print gensub("[a-z]*([0-9]*).cfg", "\\\\1", 1, $0); }'`
  echo $F
  ln -sf $FILE gastemper.cfg
  ./aspect.plt
  mv aspect.png aspect$F.png
done

qiv aspect*.png

