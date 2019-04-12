#!/bin/bash
#"./IF97" -test_all
declare -a unit_tests=("B23" "R1" "R1_PH" "R1_PS" "R2" "R2Meta" "B2bc" "R2_PH" "R2_PS" "R3" "R4" "R5")
for i in "${unit_tests[@]}"
do
    cmp -s UnitTest/$i.dat UnitTest/expected/$i.dat
    error=$?
    if [ $error -eq 0 ]
    then
       #echo "R1 and R1 are the same file"
       printf "%6b" $i "................" "\e[0;32m[OK]\e[m" "\n"
    elif [ $error -eq 1 ]
    then
       #echo "R1 and R1 differ"
       printf "%6b" $i "............" "\e[0;31m[FAILED]\e[m" "\n"
    else
       #echo "There was something wrong with the diff command"
       printf "%6b" $i "........" "\e[1;31m[DIFF ERROR]\e[m" "\n"
    fi
done
