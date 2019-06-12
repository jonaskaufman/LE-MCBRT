#!/bin/bash
#SBATCH -p gpu-shared
#SBATCH --gres=gpu:p100:1
#SBATCH --job-name="le"
#SBATCH --output="le_test_1.%j.%N.out"
#SBATCH -t 00:05:00

for N in 1000 # add elements if you want
do
    for primary_rays in 1000000
    do
        for d in C R
        do
            for M in 10 50 100 500 1000 
            do
            test_dir="le_N_${N}_M_${M}_rays_${primary_rays}_density_${d}"
            mkdir $test_dir
            cd $test_dir
            ../le-run.x $N $M $primary_rays $d 0.4 > std.out # if random then density is just ignored
            cd ..
            done
        done
    done
done

