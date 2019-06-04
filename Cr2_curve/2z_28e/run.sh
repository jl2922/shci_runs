#!/bin/bash
#SBATCH -N 16
#SBATCH -p RM
#SBATCH -c 28
##SBATCH -x xl00[1-4]
#SBATCH -t 2-00:00:00

echo $SLURM_JOB_NODELIST
echo $SLURM_CPUS_PER_TASK

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

HOSTFILE=/tmp/hosts.$SLURM_JOB_ID

srun hostname -s > $HOSTFILE

mpirun_rsh -hostfile $HOSTFILE -np 16 MV2_ENABLE_AFFINITY=0 /home/junhao/shci/shci_128

