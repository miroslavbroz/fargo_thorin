#!/bin/sh

for FILE in gaspebblevtheta*.cfg ; do
  echo $FILE
  ln -sf $FILE gaspebblevtheta.cfg
  ./gaspebblevtheta.plt
  mv gaspebblevtheta.png $FILE.png
done

qiv gaspebblevtheta*.png

