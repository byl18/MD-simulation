module load gromacs
cycle_dir=$1
cd $cycle_dir
for replica_dir in `ls -1 | grep replica | sort`;
do
    cd $replica_dir
    if [ -e prd.xtc ];
    then
    gmx_mpi trjconv -f prd.xtc -s prd.tpr -n ../../../../2_md-preparation/gleap.ndx -o prd.target.xtc <<EOF
Target
EOF
    else
    gmx_mpi trjconv -f prd_mdtraj_fit.xtc -s prd.tpr -n ../../../../2_md-preparation/gleap.ndx -o prd.target.xtc <<EOF
Target
EOF
    fi
    cd ..
done

concated=" "
for filename in `find -name prd.target.xtc | sort`;
do
    concated=" $concated $filename "
done
gmx_mpi trjcat -f $concated -cat -o prd.target.trjcat-cycle.xtc
cd ..



# gmx_mpi trjcat -f 1/cmd.skip10.whole.target.xtc 2/cmd.skip10.whole.target.xtc 3/cmd.skip10.whole.target.xtc 4/cmd.skip10.whole.target.xtc 5/cmd.skip10.whole.target.xtc -cat -o cmd.skip10.whole.target.trjcat.xtc
