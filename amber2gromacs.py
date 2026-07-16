import sys
import parmed as pmd

amber_name = sys.argv[1]
gromacs_name = sys.argv[2]

amber = pmd.load_file(f"{amber_name}.prmtop", f"{amber_name}.inpcrd")
amber.save(f"{gromacs_name}.gro")
amber.save(f"{gromacs_name}.top")
