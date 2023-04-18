#!/bin/bash                                                                                             

#SBATCH --account=NN9407K                                                                               
#SBATCH --time=1:00:00                                                                          
#SBATCH --partition=normal                                                                              
#SBATCH --nodes=16 --ntasks=128 --cpus-per-task=1                                                       

srun -n 128 ./gevolution -n 16 -m 8 -s settings.ini
