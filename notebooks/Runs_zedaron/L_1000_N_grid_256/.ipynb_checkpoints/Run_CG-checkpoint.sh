#!/bin/bash

models=("0.01" "0.1" "0.5" "1" "2")
dir="/mn/stornext/u3/hassanif/MG_evolution_project/QCG_runs_simulations_tests/L_1000_N_grid_256"

for model in "${models[@]}"
    do
    newdir="$dir/CG_ks$model"
    mkdir -p $newdir
    
    cp ./*.dat ./gevolution run.sh "$newdir"
    
    sed -e "s#k_s = xxx#k_s = $model#g" settings.ini > "$newdir"/settings.ini
    mkdir -p "$newdir/output"
    chmod 777 "$newdir/gevolution"
        # Open a new detached screen, change directory, and run the script
       # /usr/bin/screen -S "$name$name2" -dm /bin/bash -c "cd ./$name/$name2 && module load python && bash run.sh"



        # Detach the screen
#        /usr/bin/screen -S "$name$name2" -p 0 -X detach
         #2>&1 | tee -a "$log_file"  # Append both stdout and stderr to the log file

    done
