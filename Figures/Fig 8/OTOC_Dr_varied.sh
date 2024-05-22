#!/bin/bash
make Heisenberg_OTOC_Dr_varied
mkdir OTOC_driving
cd OTOC_driving
dj="0.35"
tau="0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0"
for j in $dj
do
mkdir $j 
cd $j
for i in $tau
do
mkdir $i
cp ../../Heisenberg_OTOC_Dr_varied ./$i/
cd $i
cat > tau.txt << EOF
$i
EOF
cat > dj.txt << EOF
$j
EOF
./Heisenberg_OTOC_Dr_varied
cd ..
done
cd ..
done