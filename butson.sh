#!/bin/sh
# 
# AUTHOR : ARNAUD MALLEN (mallen.arnaud@gmail.com)
# Dynamically generate butson matrix searching program
# The program has to be generated dynamically since the matrix
# searching algorithm is templated on the matrix dimensions
# 
# Prerequisites : cmake 2.6, C++11
# 
# 

printUsage()
{
	echo "butson.sh l n"
	echo "  inputs : "
	echo "    l : search matrix composed of l-th root of 1"
	echo "    n : number of L-cycles contained in each line of the matrix"
	echo "  outputs : "
	echo "    for each solution i, print i in standard ouput and write ascii file butson_l_n_i.txt containing the matrix coefficients"
}

if [ $# -ne 2 ]
then
	printUsage
	exit 1
fi

cp main.cpp ~main.cpp

OLD="LVALUE";
NEW=${1};
sed "s/$OLD/$NEW/g" "main.cpp" > "_main.cpp" && mv "_main.cpp" "main.cpp" && 
{
OLD="NVALUE";
NEW=${2};
sed "s/$OLD/$NEW/g" "main.cpp" > "_main.cpp" && mv "_main.cpp" "main.cpp"
}
# Generation
cmake .
make;

# Launch
./butson $1 $2

# Enable relaunch
mv ~main.cpp main.cpp
