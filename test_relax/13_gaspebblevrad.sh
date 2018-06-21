#!/bin/sh

for FILE in gaspebblevrad*.cfg ; do
  echo $FILE
  ln -sf $FILE gaspebblevrad.cfg
  ./gaspebblevrad.plt
  mv gaspebblevrad.png $FILE.png
done

qiv gaspebblevrad*.png

