#!/bin/bash
make Heisenberg_OTOC_Dr_varied
mkdir OTOCs
cd OTOCs
tau="0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0"
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
