#!/bin/bash
#SBATCH -N 2
#SBATCH -p LM
#SBATCH --mem=3000GB
#SBATCH -c 80
#SBATCH -x xl00[1-4]
#SBATCH -t 3-00:00:00

echo $SLURM_JOB_NODELIST
echo $SLURM_CPUS_PER_TASK

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

HOSTFILE=/tmp/hosts.$SLURM_JOB_ID

srun hostname -s > $HOSTFILE

mpirun_rsh -hostfile $HOSTFILE -np 2 MV2_ENABLE_AFFINITY=0 /home/junhao/shci/shci_128

