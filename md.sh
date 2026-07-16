######## 1 ########
mkdir 1_modeling
grep -E '^(ATOM|HETATM)' 1_3q50.pdb # delete SO4

antechamber -fi mol2 -i 2_3q50_B_PRF_opt.mol2 -fo gcrt -o prf.gau -nc 0
g16 < prf.gau > prf.out
antechamber -fi gout -i prf.out -fo ac -o prf.ac -c resp


antechamber -fi ac -i prf.ac -fo ac -o prf-renamed.ac -fa mol2 -a 2_3q50_B_PRF_opt.mol2 -ao name -at gaff2 -rn PRF
antechamber -fi ac -i prf-renamed.ac -fo mol2 -o prf.mol2
parmchk2 -f mol2 -i prf.mol2 -o prf.frcmod



######## 2 ########
mkdir 2_md-preparation
cd 2_md-preparation
pdb4amber -i ../1_modeling/2_3q50.pdb -o sys.pdb
tleap

```
source leaprc.water.opc
source leaprc.RNA.OL3
source leaprc.gaff2
loadAmberParams ../1_modeling/prf.frcmod
PRF = loadmol2 ../1_modeling/prf.mol2
saveamberparm PRF sys.inpcrd sys.prmtop
saveoff PRF ../PRF.lib

# loadAmberParams sys.prmtop
loadoff ../PRF.lib

sys = loadpdb sys.pdb
savepdb sys leap.pdb


solvatebox sys OPCBOX {23.6 19.18 15}

charge sys
addions sys K+ 53
addions sys Cl- 21
charge sys

saveamberparm sys leap.prmtop leap.inpcrd
savepdb sys leap.pdb

quit
```


python3 ../amber2gromacs.py leap gleap  # Transfer amber file into gromacs file
gmx_mpi make_ndx -f gleap.gro -o gleap.ndx # Target, NonTarget
r 1-34
name 10 Target
! 10
name 11 NonTarget
r 1-33
name 12 RNA
q

echo Target | gmx_mpi trjconv -s gleap.gro -f gleap.gro -n gleap.ndx -o gleap.target.gro
selections=( "RNA" "PRF")
for selection in ${selections[@]};
do
    echo $selection | gmx_mpi trjconv -s gleap.gro -f gleap.gro -n gleap.ndx -o $selection.gro
    gmx_mpi make_ndx -f $selection.gro -o $selection.ndx << EOF
0 & ! a H*
q
EOF
    for ((force=100; force<=1000; force+=100));
    do
        gmx_mpi genrestr -f $selection.gro -fc $force -n $selection.ndx -o posre_${selection}_${force}.itp << EOF
System
EOF
        echo "#ifdef POSRES_${force}" >> posre_${selection}.top-directive
        echo '#include "posre_'${selection}'_'$force'.itp"' >> posre_${selection}.top-directive
        echo "#endif" >> posre_${selection}.top-directive

#         gmx_mpi genrestr -f $selection.gro -fc $force -n $selection.ndx -o posre_${selection}_ca_${force}.itp << EOF
# C-alpha
# EOF
#         echo "#ifdef POSRES_${force}" >> posre_${selection}_ca.top-directive
#         echo '#include "posre_'${selection}'_ca_'$force'.itp"' >> posre_${selection}_ca.top-directive
#         echo "#endif" >> posre_${selection}_ca.top-directive

#         gmx_mpi genrestr -f $selection.gro -fc $force -n $selection.ndx -o posre_${selection}_backbone_${force}.itp << EOF
# Backbone
# EOF
#         echo "#ifdef POSRES_${force}" >> posre_${selection}_backbone.top-directive
#         echo '#include "posre_'${selection}'_backbone'$force'.itp"' >> posre_${selection}_backbone.top-directive
#         echo "#endif" >> posre_${selection}_backbone.top-directive

#         gmx_mpi genrestr -f $selection.gro -fc $force -n $selection.ndx -o posre_${selection}_heavyatoms_${force}.itp << EOF
# System_&_!H
# EOF
#         echo "#ifdef POSRES_${force}" >> posre_${selection}_heavyatoms.top-directive
#         echo '#include "posre_'${selection}'_heavyatoms_'$force'.itp"' >> posre_${selection}_heavyatoms.top-directive
#         echo "#endif" >> posre_${selection}_heavyatoms.top-directive
    done
done
# python3 $WORKSPACE_DIR/programs/insert-restraints.py 2 posre_ProteinNmr.top-directive 3 posre_RnaNmr.top-directive
gmx_mpi editconf -f gleap.gro -o gleap.pdb
cd ..


######## 3 ########
mkdir 3_min-1
cd 3_min-1
cp ../min-1.mdp .
gmx_mpi grompp -f min-1.mdp -c ../2_md-preparation/gleap.gro -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -r ../2_md-preparation/gleap.gro -o min.tpr
gmx_mpi mdrun -deffnm min
echo Potential | gmx_mpi energy -f min.edr -o Potential.xvg
gmx_mpi trjconv -s min.tpr -f min.gro -n ../2_md-preparation/gleap.ndx -o min.whole.gro -pbc whole <<EOF
System
EOF
gmx_mpi trjconv -s min.tpr -f min.whole.gro -n ../2_md-preparation/gleap.ndx -o min.whole.cluster.gro -pbc cluster <<EOF
Target
System
EOF
gmx_mpi trjconv -s min.tpr -f min.whole.cluster.gro -n ../2_md-preparation/gleap.ndx -o min.whole.cluster.mol.gro -pbc mol -center <<EOF
Target
System
EOF
gmx_mpi trjconv -s min.tpr -f min.whole.cluster.mol.gro -n ../2_md-preparation/gleap.ndx -o min.whole.cluster.mol.target.gro <<EOF
Target
EOF
gmx_mpi editconf -f min.gro -o min.pdb
cd ..


######## 4 ########
mkdir 4_nvt-1
cd 4_nvt-1
cp ../nvt-1.mdp .
gmx_mpi grompp -f nvt-1.mdp -c ../3_min-1/min.gro -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -r ../2_md-preparation/gleap.gro -o nvt.tpr
gmx_mpi mdrun -deffnm nvt
echo Temperature | gmx_mpi energy -f nvt.edr -o Temperature.xvg
gmx_mpi trjconv -s nvt.tpr -f nvt.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.gro -pbc whole <<EOF
System
EOF
gmx_mpi trjconv -s nvt.tpr -f nvt.whole.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.gro -pbc cluster <<EOF
Target
System
EOF
gmx_mpi trjconv -s nvt.tpr -f nvt.whole.cluster.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.mol.gro -pbc mol -center <<EOF
Target
System
EOF
gmx_mpi trjconv -s nvt.tpr -f nvt.whole.cluster.mol.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.mol.target.gro <<EOF
Target
EOF
gmx_mpi editconf -f nvt.gro -o nvt.pdb
cd ..


mkdir 4_nvt-2
cd 4_nvt-2
cp ../nvt-2.mdp .
gmx_mpi grompp -f nvt-2.mdp -c ../4_nvt-1/nvt.gro -t ../4_nvt-1/nvt.cpt -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -r ../2_md-preparation/gleap.gro -o nvt.tpr -maxwarn 1
gmx_mpi mdrun -deffnm nvt
echo Temperature | gmx_mpi energy -f nvt.edr -o Temperature.xvg
gmx_mpi trjconv -s nvt.tpr -f nvt.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.gro -pbc whole <<EOF
System
EOF
gmx_mpi trjconv -s nvt.tpr -f nvt.whole.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.gro -pbc cluster <<EOF
Target
System
EOF
gmx_mpi trjconv -s nvt.tpr -f nvt.whole.cluster.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.mol.gro -pbc mol -center <<EOF
Target
System
EOF
gmx_mpi trjconv -s nvt.tpr -f nvt.whole.cluster.mol.gro -n ../2_md-preparation/gleap.ndx -o nvt.whole.cluster.mol.target.gro <<EOF
Target
EOF
gmx_mpi editconf -f nvt.gro -o nvt.pdb
cd ..


nohup bash step5_npt.sh > step5_npt.log 2>&1 &

for ((i=1; i<=10; i++));
do
    if [ $i -eq 1 ];
    then
        prev_run=../4_nvt-2/nvt
    else
        prev_run=../5_npt-$(($i-1))/npt
    fi
    force=$(( (11-$i) * 100 ))
    mkdir 5_npt-$i
    cd 5_npt-$i
    cp ../npt-1.mdp .
    sed -i "s/DPOSRES_XXXX/DPOSRES_${force}/g" npt-1.mdp
    gmx_mpi grompp -f npt-1.mdp -c $prev_run.gro -t $prev_run.cpt -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -r ../2_md-preparation/gleap.gro -o npt.tpr -maxwarn 1
    gmx_mpi mdrun -deffnm npt
    echo Pressure | gmx_mpi energy -f npt.edr -o Pressure.xvg
    echo Density | gmx_mpi energy -f npt.edr -o densDensityity.xvg
    echo Temperature | gmx_mpi energy -f npt.edr -o Temperature.xvg
    gmx_mpi trjconv -s npt.tpr -f npt.gro -n ../2_md-preparation/gleap.ndx -o npt.whole.gro -pbc whole <<EOF
System
EOF
    gmx_mpi trjconv -s npt.tpr -f npt.whole.gro -n ../2_md-preparation/gleap.ndx -o npt.whole.cluster.gro -pbc cluster <<EOF
Target
System
EOF
    gmx_mpi trjconv -s npt.tpr -f npt.whole.cluster.gro -n ../2_md-preparation/gleap.ndx -o npt.whole.cluster.mol.gro -pbc mol -center <<EOF
Target
System
EOF
    gmx_mpi trjconv -s npt.tpr -f npt.whole.cluster.mol.gro -n ../2_md-preparation/gleap.ndx -o npt.whole.cluster.mol.target.gro <<EOF
Target
EOF
    gmx_mpi editconf -f npt.gro -o npt.pdb
    cd ..
done


######## 5 ########


mkdir 6_cmd-1
cd 6_cmd-1
cp ../cmd.mdp .
gmx_mpi grompp -f cmd.mdp -c ../5_npt-10/npt.gro -t ../5_npt-10/npt.cpt -p ../2_md-preparation/gleap.top -n ../2_md-preparation/gleap.ndx -o cmd.tpr -maxwarn 1
gmx_mpi dump -s cmd.tpr > cmd.1.txt
gmx_mpi mdrun -deffnm cmd
# nohup gmx_mpi mdrun -deffnm cmd -cpi cmd.cpt -append > 6.2.log 2>&1 &
# nohup gmx_mpi mdrun -deffnm cmd > 6.log 2>&1 &


echo Pressure | gmx_mpi energy -f cmd.edr -o Pressure.xvg
echo Density | gmx_mpi energy -f cmd.edr -o Density.xvg
echo Temperature | gmx_mpi energy -f cmd.edr -o Temperature.xvg

gmx_mpi trjconv -s cmd.tpr -f cmd.gro -n ../2_md-preparation-2/gleap.ndx -o cmd.whole.gro -pbc whole <<EOF
System
EOF
gmx_mpi trjconv -s cmd.tpr -f cmd.whole.gro -n ../2_md-preparation-2/gleap.ndx -o cmd.whole.target.gro <<EOF
Target
EOF

gmx_mpi trjconv -s cmd.tpr -f cmd.xtc -n ../2_md-preparation-2/gleap.ndx -o cmd.skip10.whole.xtc -skip 10 -pbc whole <<EOF
System
EOF
gmx_mpi trjconv -s cmd.tpr -f cmd.skip10.whole.xtc -n ../2_md-preparation-2/gleap.ndx -o cmd.skip10.whole.target.xtc <<EOF
Target
EOF

gmx_mpi editconf -f cmd.gro -o cmd.pdb
cd ..



### checkpoint

gmx_mpi trjconv -s cmd.tpr -f cmd.xtc -dump 100000000000000 -o tmp.gro -n ../2_md-preparation/gleap.ndx -pbc whole
gmx_mpi trjconv -s cmd.tpr -f tmp.gro -o tmp.gro -n ../2_md-preparation/gleap.ndx -pbc cluster # if tmp.gro is seperate
gmx_mpi trjconv -s cmd.tpr -f cmd.xtc -dump 100000000000000 -o tmp_f2.gro -n ../2_md-preparation/gleap.ndx -pbc whole

### center
gmx_mpi trjconv -f /work19/yilan/projects/PreQ1/6_cmd-1/cmd.xtc -n ../2_md-preparation/gleap.ndx -o center.xtc -center
gmx_mpi trjconv -f center.xtc -s cmd.tpr -n ../2_md-preparation/gleap.ndx -o nopbc.xtc -pbc mol

#### one-step
gmx_mpi trjconv -f /work19/yilan/projects/PreQ1/6_cmd-1/cmd.xtc -s /work19/yilan/projects/PreQ1/6_cmd-1/cmd.tpr -n ../2_md-preparation/gleap.ndx -o nopbc.xtc -center -pbc mol