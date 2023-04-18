#!/bin/bash

H0rc=(5 1)
k_star_values_5=(0.3 0.5 0.75 1)
k_star_values_1=(1.5 1.7 2 2.2 2.5)

for H0rc_value in "${H0rc[@]}"; do
  if [ "$H0rc_value" -eq 5 ]; then
    k_star_values=( "${k_star_values_5[@]}" )
  else
    k_star_values=( "${k_star_values_1[@]}" )
  fi

  for k_star_value in "${k_star_values[@]}"; do
    output_dir='H0rc_'$H0rc_value'_ks_'$k_star_value
    echo $output_dir;
    [ ! -d "$output_dir" ] && mkdir -p "$output_dir"
    sed -e 's/H0r_c = 1.0/H0r_c = '$H0rc_value'/g' -e 's/k_screen = 2.3/k_screen = '$k_star_value'/g' -e 's/output path = output/output path = '$output_dir'\//g' elephant.ini> settings_$output_dir.ini
    sed 's/-s settings.ini/-s settings_'$output_dir'.ini/g' cluster_run.sh > run_$output_dir.sh
    sbatch run_$output_dir.sh
  done
  echo
done
