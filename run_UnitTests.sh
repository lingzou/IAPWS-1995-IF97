#!/bin/bash
declare -a unit_tests=("B23" "R1" "R1_PH" "R1_PS" "R2" "R2Meta" "B2bc" "R2_PH" "R2_PS" "R3" "R4" "R5"
  "R3_sat_line" "R4_sat_line" "R3_rho_pT" "R3_Tx_ph" "R3_Tx_ps" "SurfTension" "Viscosity"
  "ThermCondNoEnhancement" "ThermCondR1" "ThermCondR2" "ThermCondR3" "IAPWS95_test1" "IAPWS95_test2")
ALL_CASES_SUCCEEDED=true
result=""

# remove already existing unit tests output files
for case in "${unit_tests[@]}"
do
    rm UnitTest/$case.dat
done

# run all cases
"./IF97-exe" -test_all

# do diff for each case
for case in "${unit_tests[@]}"
do
    temp=$(diff UnitTest/$case.dat UnitTest/expected/$case.dat)
    error=$?
    if [ $error -eq 0 ]
    then
       printf "%24b" $case "................" "\e[0;32m[OK]\e[m" "\n"
    elif [ $error -eq 1 ]
    then
       ALL_CASES_SUCCEEDED=false
       printf "%24b" $case "............" "\e[0;31m[FAILED]\e[m" "\n"
    else
       ALL_CASES_SUCCEEDED=false
       printf "%24b" $case "........" "\e[1;31m[DIFF ERROR]\e[m" "\n"
    fi
    result=$result$temp
done

# check diff status
if [ "$ALL_CASES_SUCCEEDED" = true ]
then
      printf "\e[0;32m[++++++++ All Unit Tests Succeed ++++++++]\e[m\n"
else
      printf "\e[0;31m[++++++ There Are Failed Unit Tests ++++++]\e[m\n"
      printf "\e[1;33m$result\e[m\n"
fi
