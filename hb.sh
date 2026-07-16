#!/bin/sh
#$ -cwd
#$ -l node_o=1
#$ -l h_rt=24:00:00
#$ -N Q1.1

module load gromacs 

for cnt in { 01 02 03 04 05 06 07 08 09 10 }; do
gmx_mpi hbond -f trial$cnt/6_cmd-1/rot_trans_fit1.xtc -s trial$cnt/6_cmd-1/cmd1.tpr -n trial$cnt/2_md-preparation/gleap.ndx -num hbxvg/hbnum$cnt.xvg <<EOF
PRF
RNA
EOF
done


