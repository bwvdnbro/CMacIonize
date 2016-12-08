#! /bin/bash

if [ $# -lt 1 ]
then
echo "Usage: ./format_script <clang-format command>"
exit
fi

files=( src/*.cpp src/*.hpp src/*.hpp.in test/*.cpp test/*.hpp python/*.cpp )

for f in "${files[@]}"
do $1 -style=file -i $f
done
