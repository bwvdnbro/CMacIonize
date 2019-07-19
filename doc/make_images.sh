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
  pdflatex -shell-escape $f
  cd ..
  filename="${f%.*}"
  mv tmp/$filename.png $filename.png
done
rm -rf tmp
