module load gromacs
# ######## 3 ########
# mkdir 3_min-1
# cd 3_min-1
# cp ../min-1.mdp .
# gmx_mpi grompp -f min-1.mdp -c ../2_md-preparation/gleap.gro -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -r ../2_md-preparation/gleap.gro -o min.tpr
# gmx_mpi mdrun -deffnm min
# echo Potential | gmx_mpi energy -f min.edr -o Potential.xvg
# gmx_mpi trjconv -s min.tpr -f min.gro -n ../2_md-preparation/gleap.ndx -o min.whole.gro -pbc whole <<EOF
# System
# EOF
# gmx_mpi trjconv -s min.tpr -f min.whole.gro -n ../2_md-preparation/gleap.ndx -o min.whole.cluster.gro -pbc cluster <<EOF
# Target
# System
# EOF
# gmx_mpi trjconv -s min.tpr -f min.whole.cluster.gro -n ../2_md-preparation/gleap.ndx -o min.whole.cluster.mol.gro -pbc mol -center <<EOF
# Target
# System
# EOF
# gmx_mpi trjconv -s min.tpr -f min.whole.cluster.mol.gro -n ../2_md-preparation/gleap.ndx -o min.whole.cluster.mol.target.gro <<EOF
# Target
# EOF
# gmx_mpi editconf -f min.gro -o min.pdb
# cd ..


# ######## 4 ########
# mkdir 4_nvt-1
# cd 4_nvt-1
# cp ../nvt-1.mdp .
# gmx_mpi grompp -f nvt-1.mdp -c ../3_min-1/min.gro -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -r ../2_md-preparation/gleap.gro -o nvt.tpr
# gmx_mpi mdrun -deffnm nvt
# echo Temperature | gmx_mpi energy -f nvt.edr -o Temperature.xvg
# gmx_mpi trjconv -s nvt.tpr -f nvt.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.gro -pbc whole <<EOF
# System
# EOF
# gmx_mpi trjconv -s nvt.tpr -f nvt.whole.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.gro -pbc cluster <<EOF
# Target
# System
# EOF
# gmx_mpi trjconv -s nvt.tpr -f nvt.whole.cluster.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.mol.gro -pbc mol -center <<EOF
# Target
# System
# EOF
# gmx_mpi trjconv -s nvt.tpr -f nvt.whole.cluster.mol.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.mol.target.gro <<EOF
# Target
# EOF
# gmx_mpi editconf -f nvt.gro -o nvt.pdb
# cd ..


# mkdir 4_nvt-2
# cd 4_nvt-2
# cp ../nvt-2.mdp .
# gmx_mpi grompp -f nvt-2.mdp -c ../4_nvt-1/nvt.gro -t ../4_nvt-1/nvt.cpt -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -r ../2_md-preparation/gleap.gro -o nvt.tpr -maxwarn 1
# gmx_mpi mdrun -deffnm nvt
# echo Temperature | gmx_mpi energy -f nvt.edr -o Temperature.xvg
# gmx_mpi trjconv -s nvt.tpr -f nvt.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.gro -pbc whole <<EOF
# System
# EOF
# gmx_mpi trjconv -s nvt.tpr -f nvt.whole.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.gro -pbc cluster <<EOF
# Target
# System
# EOF
# gmx_mpi trjconv -s nvt.tpr -f nvt.whole.cluster.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.mol.gro -pbc mol -center <<EOF
# Target
# System
# EOF
# gmx_mpi trjconv -s nvt.tpr -f nvt.whole.cluster.mol.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.mol.target.gro <<EOF
# Target
# EOF
# gmx_mpi editconf -f nvt.gro -o nvt.pdb
# cd ..



# for ((i=1; i<=10; i++));
# do
#     if [ $i -eq 1 ];
#     then
#         prev_run=../4_nvt-2/nvt
#     else
#         prev_run=../5_npt-$(($i-1))/npt
#     fi
#     force=$(( (11-$i) * 100 ))
#     mkdir 5_npt-$i
#     cd 5_npt-$i
#     cp ../npt-1.mdp .
#     sed -i "s/DPOSRES_XXXX/DPOSRES_${force}/g" npt-1.mdp
#     gmx_mpi grompp -f npt-1.mdp -c $prev_run.gro -t $prev_run.cpt -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -r ../2_md-preparation/gleap.gro -o npt.tpr -maxwarn 1
#     gmx_mpi mdrun -deffnm npt
#     echo Pressure | gmx_mpi energy -f npt.edr -o Pressure.xvg
#     echo Density | gmx_mpi energy -f npt.edr -o densDensityity.xvg
#     echo Temperature | gmx_mpi energy -f npt.edr -o Temperature.xvg
#     gmx_mpi trjconv -s npt.tpr -f npt.gro -n ../2_md-preparation/gleap.ndx -o npt.whole.gro -pbc whole <<EOF
# System
# EOF
#     gmx_mpi trjconv -s npt.tpr -f npt.whole.gro -n ../2_md-preparation/gleap.ndx -o npt.whole.cluster.gro -pbc cluster <<EOF
# Target
# System
# EOF
#     gmx_mpi trjconv -s npt.tpr -f npt.whole.cluster.gro -n ../2_md-preparation/gleap.ndx -o npt.whole.cluster.mol.gro -pbc mol -center <<EOF
# Target
# System
# EOF
#     gmx_mpi trjconv -s npt.tpr -f npt.whole.cluster.mol.gro -n ../2_md-preparation/gleap.ndx -o npt.whole.cluster.mol.target.gro <<EOF
# Target
# EOF
#     gmx_mpi editconf -f npt.gro -o npt.pdb
#     cd ..
# done


# ######## 5 ########


# mkdir 6_cmd-1
cd 6_cmd-1
# cp ../cmd.mdp .
# gmx_mpi grompp -f cmd.mdp -c ../5_npt-10/npt.gro -t ../5_npt-10/npt.cpt -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -o cmd.tpr -maxwarn 1
# gmx_mpi dump -s cmd.tpr > cmd.1.txt

for cnt in {1..5}; do
if [ ! -f cmd$cnt.tpr ]; then
	if [ $cnt -eq 1 ]; then
		gmx_mpi grompp -f cmd.mdp -c ../5_npt-10/npt.gro -t ../5_npt-10/npt.cpt -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -o cmd$cnt.tpr -maxwarn 1
	else
		gmx_mpi grompp -f cmd.mdp -c cmd$(($cnt-1)).gro -t cmd$(($cnt-1)).cpt -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -o cmd$cnt.tpr -maxwarn 1
	fi
fi
if [ ! -f cmd$cnt.gro ]; then
	if [ -f cmd$cnt.cpt ]; then 
		gmx_mpi mdrun -deffnm cmd$cnt -cpi cmd$cnt.cpt
	else
		gmx_mpi mdrun -deffnm cmd$cnt 
	fi
fi
echo Pressure | gmx_mpi energy -f cmd$cnt.edr -o Pressure$cnt.xvg
echo Density | gmx_mpi energy -f cmd$cnt.edr -o Density$cnt.xvg
echo Temperature | gmx_mpi energy -f cmd$cnt.edr -o Temperature$cnt.xvg


gmx_mpi trjconv -s cmd$cnt.tpr -f cmd$cnt.gro -n ../2_md-preparation/gleap.ndx -o cmd$cnt.whole.gro -pbc whole <<EOF
System
EOF
gmx_mpi trjconv -s cmd$cnt.tpr -f cmd$cnt.whole.gro -n ../2_md-preparation/gleap.ndx -o cmd$cnt.whole.target.gro <<EOF
Target
EOF

gmx_mpi trjconv -s cmd$cnt.tpr -f cmd$cnt.xtc -n ../2_md-preparation/gleap.ndx -o cmd$cnt.skip10.whole.xtc -skip 10 -pbc whole <<EOF
System
EOF
gmx_mpi trjconv -s cmd$cnt.tpr -f cmd$cnt.skip10.whole.xtc -n ../2_md-preparation/gleap.ndx -o cmd$cnt.skip10.whole.target.xtc <<EOF
Target
EOF

gmx_mpi editconf -f cmd$cnt.gro -o cmd$cnt.pdb
gmx_mpi trjconv -f cmd$cnt.xtc -s cmd$cnt.tpr -n ../2_md-preparation/gleap.ndx -o pbc$cnt.xtc -center -pbc mol <<EOF
Center
Target
EOF

gmx_mpi trjconv -f pbc$cnt.xtc -s cmd$cnt.tpr -o rot_trans_fit$cnt.xtc -fit rot+trans -n ../2_md-preparation/gleap.ndx <<EOF
Target
Target
EOF


done
cd ..
