#!/bin/sh

for FILE in gasvrad*.cfg ; do
  echo $FILE
  ln -sf $FILE gasvrad.cfg
  ./gasvrad.plt
  mv gasvrad.png $FILE.png
done

qiv gasvrad*.png

