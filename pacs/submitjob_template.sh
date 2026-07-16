#!/bin/sh
#$ -cwd
#$ -l node_q=1
#$ -l h_rt=24:00:00
#$ -N Plog_0.9
##$ -hold_jid 6573865

module purge
module load cuda/12.3.2
module load openmpi/5.0.2-gcc
module load gromacs/2024.2-plumed

# qsub -v TRIAL_ID=11 这样传进来
if [ -z "$TRIAL_ID" ]; then
    echo "Error: TRIAL_ID is not set"
    exit 1
fi

echo "Running in $(pwd)"
echo "TRIAL_ID=$TRIAL_ID"

/gs/bs/tga-KitaoLab-NexusYtec/yilan/softwares/miniconda3/envs/pacsmd/bin/pacs mdrun -t "$TRIAL_ID" -f input_gpu.toml
