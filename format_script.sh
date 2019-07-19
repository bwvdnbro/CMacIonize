#! /bin/bash

command -v clang-format-3.8 >/dev/null 2>&1 || \
  { echo >&2 "This script requires clang-format-3.8, but it is not installed!" \
             "Aborting."; exit 1; }

files=( src/*.cpp src/*.hpp src/*.cpp.in src/*.hpp.in test/*.cpp test/*.hpp \
        test/*.c python/*.cpp timing/*.cpp timing/*.hpp c/*.h )

for f in "${files[@]}"
do clang-format-3.8 -style=file -i $f
done
