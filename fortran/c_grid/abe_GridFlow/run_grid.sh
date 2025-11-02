#!/bin/sh
#PBS -l nodes=yato215:ppn=4
#PBS -N NACA_Grid
#PBS -o /dev/null
#PBS -j oe

OMP_NUM_THREADS=4; export OMP_NUM_THREADS
ulimit -s unlimited
cd $PBS_O_WORKDIR
(echo 8; echo 2) | (time ./temp/GridFlow/Grid.exe) &> Log_Grid.txt
