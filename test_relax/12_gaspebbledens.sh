#!/bin/sh

for FILE in gaspebbledens[0-9]*.cfg ; do
  F=`echo $FILE | awk '{ print gensub("[a-z]*([0-9]*).cfg", "\\\\1", 1, $0); }'`
  echo $FILE
  ln -sf $FILE gaspebbledens.cfg
  ln -sf gasdens$F.cfg gasdens.cfg
  ./gaspebbledens.plt
  mv gaspebbledens.png $FILE.png
done

qiv gaspebbledens*.png

