#!/bin/sh

for FILE in gastemper[0-9]*.cfg ; do
  echo $FILE
  ln -sf $FILE gastemper.cfg
  ./gastemper.plt
  mv gastemper.png $FILE.png
done

qiv gastemper*.png

