#!/bin/bash

module load Intel_parallel_studio/2020

mpirun -np 16 ./gevolution -n 4 -m 4 -s settings.ini > out.txt
