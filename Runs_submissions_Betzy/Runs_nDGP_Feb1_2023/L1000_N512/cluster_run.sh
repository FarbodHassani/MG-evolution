#!/bin/bash                                                                                             

#SBATCH --account=NN9407K                                                                               
#SBATCH --time=2:00:00                                                                          
#SBATCH --partition=normal                                                                              
#SBATCH --nodes=64 --ntasks=512 --cpus-per-task=1                                                       

srun -n 512 ./gevolution -n 32 -m 16 -s settings.ini
