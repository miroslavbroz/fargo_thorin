#!/bin/sh

for FILE in gasdens[0-9]*.cfg ; do
  echo $FILE
  ln -sf $FILE gasdens.cfg
  ./gasdens.plt
  mv gasdens.png $FILE.png
done

qiv gasdens*.png

