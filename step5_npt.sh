

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


