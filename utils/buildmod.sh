#!/usr/bin/bash

move_files () {

    files=$(ls -p | grep -v /)

    for file in $files
    do
        ext=${file##*.}
        
        if [ $ext = "o" ]
        then
            mv $file "build/linux/$file"
        elif [ $ext = "mod" ]
        then
            mv $file "modules/$file"
        fi
    done

}

for arg in $@ 
do
    gfortran -c src/$arg.f90
    move_files
done
