#!/bin/bash
make Heisenberg_Driven_Varied
cd run
cd dj0.5
tau="0.1 0.2 0.3 0.4"
for i in $tau
do
mkdir $i
cp ../../Heisenberg_Driven_Varied ./$i/
cd $i
cat > tau.txt << EOF
$i
EOF
./Heisenberg_Driven_Varied
cd ..
done
