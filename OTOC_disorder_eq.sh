#!/bin/bash
make Heisenberg_OTOC_disorder
mkdir OTOC_disorder
cd OTOC_disorder
dj="0.05"
for i in $dj
do
mkdir $i
cp ../Heisenberg_OTOC_disorder ./$i/
cd $i
cat > dj.txt << EOF
$i
EOF
./Heisenberg_OTOC_disorder
cd ..
done
