#!/bin/bash
make Heisenberg_OTOC_Dr_varied
mkdir OTOCsdj0.35
cd OTOCsdj0.35
tau="0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0"
for i in $tau
do
mkdir $i
cp ../Heisenberg_OTOC_Dr_varied ./$i/
cd $i
cat > tau.txt << EOF
$i
EOF
./Heisenberg_OTOC_Dr_varied
cd ..
done
