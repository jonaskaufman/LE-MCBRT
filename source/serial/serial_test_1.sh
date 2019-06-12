#!/bin/bash
#SBATCH --job-name="serial"
#SBATCH --output="serial_test_1.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH -t 00:10:00

N=1000
for rays in 10000 100000 1000000
do
    test_dir=serial_N_${N}_rays_${rays}_density_R
    mkdir $test_dir
    cd $test_dir
    ../serial-run $N $rays R > std.out
    cd ..
done
