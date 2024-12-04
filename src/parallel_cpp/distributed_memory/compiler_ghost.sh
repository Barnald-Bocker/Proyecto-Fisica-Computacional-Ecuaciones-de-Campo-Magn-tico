#!/bin/bash

g++ -fopenmp Capacitor_ghost.cpp
for i in {1..8}
do
  export OMP_NUM_THREADS=$i
  ./a.out | head -n 2
done
#./a.out | tail -n 100 > capacitor.csv
#python3 plot.py
