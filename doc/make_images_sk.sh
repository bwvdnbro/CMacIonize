#! /bin/bash

if [ $# -gt 0 ]
then
  sketch=$1
else
  echo "Usage: ./make_images_sk.sh SKETCH_PATH [FILE_PATTERN]"
  exit 1
fi

if [ $# -gt 1 ]
then
  files=( $2 )
else
  files=( *.tex )
fi

mkdir tmp
for f in "${files[@]}"
do
  filename="${f%.*}"
  cp $f tmp/$f
  cp "$filename"_sk.sk tmp/"$filename"_sk.sk
  cd tmp
  $sketch -o "$filename"_sk.tex "$filename"_sk.sk
  pdflatex -shell-escape $f
  cd ..
  mv tmp/$filename.png $filename.png
done
rm -rf tmp
