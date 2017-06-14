#! /bin/bash

if [ $# -gt 0 ]
then
  files=( $1 )
else
  files=( *.tex )
fi

mkdir tmp
for f in "${files[@]}"
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
