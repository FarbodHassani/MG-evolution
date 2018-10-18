#! /bin/bash

make
LD_LIBRARY_PATH=/home/farbod/packages/hdf5-1.8.17/hdf5/lib
export LD_LIBRARY_PATH
export OMPI_MCA_rmaps_base_oversubscribe=1
mpirun -np 8 ./gevolution -n 4 -m 2 -s settings.ini

rm gevolution

