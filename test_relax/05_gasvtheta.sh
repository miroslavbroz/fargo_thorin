#!/bin/sh

for FILE in gasvtheta*.cfg ; do
  echo $FILE
  ln -sf $FILE gasvtheta.cfg
  ./gasvtheta.plt
  mv gasvtheta.png $FILE.png
done

qiv gasvtheta*.png

