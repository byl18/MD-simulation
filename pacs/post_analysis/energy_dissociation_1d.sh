#!/bin/sh
#$ -cwd
#$ -l node_o=1
#$ -l h_rt=24:00:00
#$ -N post_pacs
##$ -hold_jid 6503659

export OPENBLAS_NUM_THREADS=8
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8

cd /gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/pacs2/post_analysis
/gs/bs/tga-KitaoLab-NexusYtec/yilan/softwares/miniconda3/envs/pacsmd/bin/python -u /gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/pacs2/post_analysis/energy_dissociation_1d_plot.py
