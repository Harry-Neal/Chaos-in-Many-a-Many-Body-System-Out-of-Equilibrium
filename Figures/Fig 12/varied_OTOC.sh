#!/bin/bash
make Heisenberg_NNN_OTOC_var
dj="0.05 0.1 0.25 0.5"
tau="0.1 0.2 0.3 0.4"
mkdir OTOC_dr
cd OTOC_dr
for j in $dj
do
mkdir $j
cd $j
for i in $tau
do
mkdir $i
cp ../../Heisenberg_NNN_OTOC_var ./$i/
cd $i
cat > tau.txt << EOF
$i
EOF
cat > dj.txt << EOF
$j
EOF
./Heisenberg_NNN_OTOC_var
cd ..
done
cd ..
done