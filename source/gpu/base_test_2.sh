#!/bin/bash
#SBATCH -p gpu-shared
#SBATCH --gres=gpu:p100:1
#SBATCH --job-name="base"
#SBATCH --output="base_test_2.%j.%N.out"
#SBATCH -t 00:10:00

for N in 5000 # add elements if you want
do
    for primary_rays in 1000000 10000000
    do
        for d in C R
        do
            test_dir="base_N_${N}_rays_${primary_rays}_density_${d}"
            mkdir $test_dir
            cd $test_dir
            ../base-run.x $N $primary_rays $d 0.4 > std.out # if random then density is just ignored
            cd ..
        done
    done
done

