#!/bin/sh
#$ -cwd
#$ -l cpu_80=1
#$ -l h_rt=24:00:00
#$ -N movie_Pacs_1.2
#$ -hold_jid 6503656

module load gromacs

/gs/bs/tga-Kitao-Lab/yilan/softwares/miniconda3/envs/pacsmd2/bin/pacs mdrun -t 2 -f input.toml



# gmx_mpi trjconv -s ../6_cmd-1/cmd5.tpr             -n ../2_md-preparation/gleap.ndx             -f ../6_cmd-1/cmd5.gro             -o ../6_cmd-1/cmd5_fit.gro             -center
# gmx_mpi trjconv -s ../6_cmd-1/cmd5.tpr             -n ../2_md-preparation/gleap.ndx             -f ../6_cmd-1/cmd5_fit.gro             -o ../6_cmd-1/cmd5_fit_pbcatom.gro             -pbc atom
