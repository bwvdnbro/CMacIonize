#! /bin/bash

mkdir tmp
for f in *.tex
do
  cp $f tmp/$f
  cd tmp
  filename="${f%.*}"
  latex $filename
  convert $filename.dvi $filename.png
  cd ..
  mv tmp/$filename.png $filename.png
done
rm -rf tmp
