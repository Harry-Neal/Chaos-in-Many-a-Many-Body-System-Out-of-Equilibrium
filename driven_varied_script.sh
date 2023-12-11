#!/bin/bash
make Heisenberg_Driven_Varied
mkdir run
cd run
tau="0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0"
for i in $tau
do
mkdir $i
cp ../Heisenberg_Driven_Varied ./$i/
cd $i
cat > tau.txt << EOF
$i
EOF
./Heisenberg_Driven_Varied
cd ..
done
