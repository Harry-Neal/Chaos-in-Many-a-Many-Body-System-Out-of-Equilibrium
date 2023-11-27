#!/bin/bash

make Heisenberg_Driven_Varied

mkdir run

cd run

tau="0.55 0.65 0.75 0.85 0.95 1.05 1.15 1.25 1.35 1.45 1.55 1.66 1.75 1.85 1.95 2.05 2.1 2.15 2.2 2.25"

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
