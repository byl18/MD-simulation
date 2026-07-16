#!/bin/sh
#$ -cwd
#$ -l node_q=1
#$ -l h_rt=24:00:00
#$ -N post_pacs
##$ -hold_jid 6503659
WORKSPACE_DIR=/gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/pacs2
module load gromacs

all_cycles=""
for trial in $(seq 11 20);
do
    trial_dir=`printf "trial%03d" "$trial"`
    cd $trial_dir

    for cycle_dir in `ls -1 | grep cycle`;
    do
        all_cycles="$all_cycles `pwd`/$cycle_dir"
    done
    cd ..
done
echo $all_cycles | sort | xargs -d " " -L 1 -P 16 sh $WORKSPACE_DIR/post_analysis/pacs-postscript-helper.sh



for trial in $(seq 11 20);
do
    (
    trial_dir=`printf "trial%03d" "$trial"`
    cd $trial_dir

    concated=" "
    for filename in `find -name prd.target.trjcat-cycle.xtc | sort`;
    do
        concated=" $concated $filename "
    done

    rm prd.target.trjcat-all.xtc
    gmx_mpi trjcat -f $concated -cat -o prd.target.trjcat-all.xtc

    find -maxdepth 2 -name prd.target.trjcat-cycle.xtc | xargs rm

    find -maxdepth 3 -name input.gro    | xargs -P1024 rm
    find -maxdepth 3 -name prd.gro      | xargs -P1024 rm
    find -maxdepth 3 -name prd.cpt      | xargs -P1024 rm
    find -maxdepth 3 -name prd_prev.cpt | xargs -P1024 rm


    cd ..
    )&
done
wait
