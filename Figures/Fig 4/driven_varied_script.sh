#!/bin/bash
make Heisenberg_Driven_Varied
mkdir Energy_absorbtion
cd Energy_absorbtion
dj="0.5"
tau="2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0"
for j in $dj
do 
mkdir $j
cd $j
for i in $tau
do
mkdir $i
cp ../../Heisenberg_Driven_Varied ./$i/
cd $i
cat > tau.txt << EOF
$i
EOF
cat > dj.txt << EOF
$j
EOF
./Heisenberg_Driven_Varied
cd ..
done 
cd ..
done
