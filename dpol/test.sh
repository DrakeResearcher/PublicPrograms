#!/usr/bin/env bash

while getopts c o; do
    case $o in
        (c) COMP="yes";;
    esac
done
shift "$((OPTIND - 1))"

MAIN="src/dbigdq.f90"
MOD=("build/linux/dqmodule.o" "build/linux/dqfun.o" "build/linux/wavext.o")
OUT="dbigdq_test"
FLAGS=("-std=legacy" "-mcmodel=medium")
INC=("modules")

TEST_IGNORE=(45)

gfortran $MAIN ${MOD[@]} -o $OUT ${FLAGS[@]} -I ${INC[@]}

if [ $? -eq 0 ]
then
    echo "copying test .dat file..."
    cp data/test.dat data/dall2020.dat
    clear

    if [ $COMP == "yes" ]; then
        rm ./dbigdq.out
        ./$OUT | tee ./dbigdq.out
        util compare.py ./dbigdq.out ./baseline.out ${TEST_IGNORE[@]}
    else
        ./$OUT
    fi
 
    
else
    echo "Compilation failed"
fi
