#!/bin/bash

g++ Capacitor.cpp
./a.out > capacitor.csv
python3 plot.py
