#!/bin/bash

mpicxx Capacitor.cpp

#mpirun -np 2 a.out | head -n 2

mpirun -np 2 a.out | tail -n +3 > capacitor.csv
python3 plot.py
