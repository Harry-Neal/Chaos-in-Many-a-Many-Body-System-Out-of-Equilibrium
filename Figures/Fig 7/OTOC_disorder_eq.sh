#!/bin/bash
make Heisenberg_OTOC_disorder
beta="1.06"
dj="0.0 0.05 0.1 0.15 0.2"
mkdir OTOC_dj_varied
cd OTOC_dj_varied
for j in $beta
do
mkdir $j
cd $j
for i in $dj
do
mkdir $i
cp ../../Heisenberg_OTOC_disorder ./$i/
cd $i
cat > beta.txt << EOF
$j
EOF
cat > dj.txt << EOF
$i
EOF
./Heisenberg_OTOC_disorder
cd ..
done
cd ..
done
