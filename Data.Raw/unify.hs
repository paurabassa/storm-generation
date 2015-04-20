#!/bin/bash

# script to unify severel files on water level and surge into a single one

for file in HistoricData/r052*.dat HistoricData/r100*.dat 
do 
tail -n +2 $file >> resid.dat
done

for file in HistoricData/s052*.dat HistoricData/s100*.dat     
do
tail -n +2 $file >> total.dat
done

