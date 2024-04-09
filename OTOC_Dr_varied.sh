#!/bin/bash
make Heisenberg_OTOC_Dr_varied
cd OTOCs
tau="0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.5 5.75 6.0"
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
