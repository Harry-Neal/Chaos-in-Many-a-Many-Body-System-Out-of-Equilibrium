#!/bin/bash

tau="0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2.0 2.05 2.1 2.25 2.2 2.25"

cd run

rm ./results.txt

for i in $tau
do
cd $i
cat Avg_energy.dat >> ../results.txt
echo >> ../results.txt
cd ..
done
