#!/bin/bash
make Heisenberg_NNN_Dr_var
mkdir energy_absorbtion
cd energy_absorbtion
tau="0.1 0.2 0.3 0.4"
for i in $tau
do
mkdir $i
cp ../Heisenberg_NNN_Dr_var ./$i/
cd $i
cat > tau.txt << EOF
$i
EOF
./Heisenberg_NNN_Dr_var
cd ..
done