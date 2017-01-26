#!/bin/bash

for DAMP in `seq 0.76 0.01 0.79`; do
    echo "running for damping of $DAMP"
    #lmp_ompi_g++ -var DAMP ${DAMP} -in run.lmp &> damp_${DAMP}.out
    mpirun -np 4 lmp_ompi_g++ -var DAMP ${DAMP} -in run.lmp &> damp_${DAMP}.out
    awk 'NF == 5' dip.${DAMP}.out > dip.${DAMP}.dat
done
