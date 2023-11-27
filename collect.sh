#!/bin/bash

tau="0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0"

cd run

rm ./results.txt

for i in $tau
do
cd $i
cat Avg_energy.dat >> ../results.txt
echo >> ../results.txt
cd ..
done
