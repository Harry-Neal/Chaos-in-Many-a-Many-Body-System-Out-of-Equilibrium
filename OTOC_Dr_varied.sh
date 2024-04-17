#!/bin/bash
make Heisenberg_OTOC_Dr_varied
mkdir OTOCsdj0.5
cd OTOCsdj0.5
tau=" 1 1.5 2 2.5 3 3.5 4"
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
