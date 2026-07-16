#!/bin/sh
#$ -cwd
#$ -l node_q=1
#$ -l h_rt=24:00:00
#$ -N Q1.1
#$ -hold_jid 5594573   

module load gromacs 

bash md_total.sh 




