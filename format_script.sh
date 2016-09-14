#! /bin/bash

if [ $# -lt 1 ]
then
echo "Usage: ./format_script <clang-format command>"
exit
fi

files=( src/*.cpp src/*.hpp test/*.cpp test/*.hpp )

for f in "${files[@]}"
do $1 -style=file -i $f
done
