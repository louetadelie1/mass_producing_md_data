#!/bin/bash -l
#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J protein
#SBATCH --time=00:30:00
#SBATCH --nodes=8
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=23

module purge
module load gcc/13 cmake/3.30 openmpi/5.0
module load gromacs/2025.1

plumed_prefix=/u/adlouet/src/plumed2
export PLUMED_KERNEL=$plumed_prefix/lib/libplumedKernel.so
export LD_LIBRARY_PATH=$plumed_prefix/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export GMX_FORCE_UPDATE_DEFAULT_GPU=1
export GMX_DISABLE_GPU_TIMING=1

srun --mpi=pmix gmx_mpi mdrun -v -ntomp $OMP_NUM_THREADS -s prod_final.tpr -multidir rep_*/prod/ -plumed ../../plumed_files/plumed.dat -maxh 00.25 -nb gpu -noappend -cpi state

