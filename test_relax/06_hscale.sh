#!/bin/sh

for FILE in gastemper*.cfg ; do
  F=`echo $FILE | awk '{ print gensub("[a-z]*([0-9]*).cfg", "\\\\1", 1, $0); }'`
  echo $F
  ln -sf $FILE gastemper.cfg
  ./hscale.plt
  mv hscale.png hscale$F.png
done

qiv hscale*.png

