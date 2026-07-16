# Molecular Dynamics (MD) Simulation Workflow for Protein-RNA-Ligand Systems

## Overview

This document describes a complete MD simulation pipeline for studying protein-RNA-ligand complexes, specifically designed for systems like the PreQ1 riboswitch with bound ligands (SAM, SAH, EEM). The workflow covers everything from ligand parameterization to production MD simulations and trajectory processing.



## 1. Ligand Parameterization (Step 1)

### Purpose

Prepare force field parameters for non-standard ligands (small molecules) that are not available in standard force fields.

### Process

```bash
# Remove sulfate ions from crystal structure
grep -E '^(ATOM|HETATM)' 1_3q50.pdb

# Generate Gaussian input for RESP charge calculation
antechamber -fi mol2 -i ligand_opt.mol2 -fo gcrt -o ligand.gau -nc 0

# Run Gaussian to calculate electrostatic potential
g16 < ligand.gau > ligand.out

# Extract RESP charges
antechamber -fi gout -i ligand.out -fo ac -o ligand.ac -c resp

# Rename atoms to match original structure and assign GAFF2 parameters
antechamber -fi ac -i ligand.ac -fo ac -o ligand-renamed.ac -fa mol2 \
            -a ligand_opt.mol2 -ao name -at gaff2 -rn LIG

# Convert to mol2 format
antechamber -fi ac -i ligand-renamed.ac -fo mol2 -o ligand.mol2

# Generate force field modification file
parmchk2 -f mol2 -i ligand.mol2 -o ligand.frcmod
```

### Key Outputs

- `ligand.mol2`: Ligand structure with charges
- `ligand.frcmod`: Force field corrections for the ligand



## 2. System Preparation (Step 2)

### Purpose

Build the complete simulation system including protein, RNA, ligand, ions, and water.

### Key Steps

#### A. Prepare Protein Structure

```bash
pdb4amber -i ../1_modeling/protein.pdb -o sys.pdb
```

This removes water molecules, adds missing atoms, and standardizes the PDB file.

#### B. Load Parameters in LEaP

```tleap
source leaprc.water.opc          # OPC water model
source leaprc.RNA.OL3            # RNA parameters
source leaprc.gaff2              # General Amber force field
loadAmberParams ligand.frcmod    # Custom ligand parameters
PRF = loadmol2 ligand.mol2       # Load ligand
saveamberparm PRF sys.inpcrd sys.prmtop
saveoff PRF ../PRF.lib           # Save as library for future use
```

#### C. Build the Complex

```tleap
loadoff ../PRF.lib
sys = loadpdb sys.pdb
solvatebox sys OPCBOX {23.6 19.18 15}  # Solvate with water box
charge sys
addions sys K+ 53    # Add counterions
addions sys Cl- 21   # Neutralize system
saveamberparm sys leap.prmtop leap.inpcrd
```

#### D. Convert to GROMACS Format

```bash
python3 ../amber2gromacs.py leap gleap  # Convert Amber to GROMACS
```

#### E. Generate Index Files

```bash
gmx_mpi make_ndx -f gleap.gro -o gleap.ndx
# Define groups: Target (ligand), NonTarget (protein+RNA), RNA
# Example: r 1-34 → Target, r 1-33 → RNA
```

#### F. Position Restraints for Equilibration

For each component (RNA and ligand), generate position restraints at different force constants:

```bash
gmx_mpi genrestr -f selection.gro -fc 100 -n selection.ndx -o posre_selection_100.itp
# Forces: 100, 200, 300, ..., 1000 kJ/mol/nm²
```



## 3. Energy Minimization (Step 3)

### Purpose

Remove bad contacts and relax the initial structure.

### Process

```bash
gmx_mpi grompp -f min-1.mdp -c gleap.gro -p gleap.top -n gleap.ndx -r gleap.gro -o min.tpr
gmx_mpi mdrun -deffnm min
```

### Check Energy

```bash
echo Potential | gmx_mpi energy -f min.edr -o Potential.xvg
```

### Post-processing

Remove PBC artifacts and center the target molecule:

```bash
gmx_mpi trjconv -s min.tpr -f min.gro -n gleap.ndx -o min.whole.gro -pbc whole
gmx_mpi trjconv -s min.tpr -f min.whole.gro -n gleap.ndx -o min.whole.cluster.gro -pbc cluster
gmx_mpi trjconv -s min.tpr -f min.whole.cluster.gro -n gleap.ndx -o min.whole.cluster.mol.gro -pbc mol -center
```



## 4. Equilibration (Step 4)

### Purpose

Gradually heat the system to target temperature (300 K) and equilibrate density.

### NVT Ensemble (Constant Volume)

Two stages:

1. **NVT-1**: Initial heating to 300 K
2. **NVT-2**: Continued equilibration at constant volume

```bash
# NVT-1
gmx_mpi grompp -f nvt-1.mdp -c ../3_min-1/min.gro -p gleap.top -n gleap.ndx -r gleap.gro -o nvt.tpr
gmx_mpi mdrun -deffnm nvt

# NVT-2
gmx_mpi grompp -f nvt-2.mdp -c ../4_nvt-1/nvt.gro -t ../4_nvt-1/nvt.cpt -p gleap.top -n gleap.ndx -r gleap.gro -o nvt.tpr
gmx_mpi mdrun -deffnm nvt
```

### Check Temperature

```bash
echo Temperature | gmx_mpi energy -f nvt.edr -o Temperature.xvg
```

### NPT Ensemble (Constant Pressure)

Gradually reduce position restraints over 10 steps (force constants from 1000 to 100 kJ/mol/nm²):

```bash
for ((i=1; i<=10; i++)); do
    force=$(( (11-$i) * 100 ))  # 1000, 900, ..., 100
    sed -i "s/DPOSRES_XXXX/DPOSRES_${force}/g" npt-1.mdp
    gmx_mpi grompp -f npt-1.mdp -c prev_run.gro -t prev_run.cpt -p gleap.top -n gleap.ndx -o npt.tpr
    gmx_mpi mdrun -deffnm npt
done
```



## 5. Production MD (Step 5)

### Purpose

Run long MD simulations to sample conformational space and collect data for binding free energy calculations.

### Multiple Independent Runs

Use checkpointing to run long simulations:

```bash
# Grompp
gmx_mpi grompp -f cmd.mdp -c ../5_npt-10/npt.gro -t ../5_npt-10/npt.cpt -p gleap.top -n gleap.ndx -o cmd1.tpr

# Run MD (can restart from checkpoint)
gmx_mpi mdrun -deffnm cmd1  # First run
gmx_mpi mdrun -deffnm cmd1 -cpi cmd1.cpt  # Continue from checkpoint
```

### Multiple Replicas

Five independent production runs:

```bash
for cnt in {1..5}; do
    # Grompp and mdrun for each replica
    # Each run continues from previous checkpoint
done
```



## 6. Trajectory Processing

### Purpose

Prepare trajectories for binding free energy calculations (PCA, clustering, MSM).

### Remove PBC Artifacts and Center

```bash
gmx_mpi trjconv -f cmd.xtc -s cmd.tpr -n gleap.ndx -o nopbc.xtc -center -pbc mol
```

### Extract Target (Ligand) Trajectory

```bash
gmx_mpi trjconv -s cmd.tpr -f cmd.xtc -n gleap.ndx -o target.xtc -pbc mol
```

### Fit to Target for Rotation/Translation Removal

```bash
gmx_mpi trjconv -f nopbc.xtc -s cmd.tpr -o rot_trans_fit.xtc -fit rot+trans -n gleap.ndx
```

### Downsample Trajectories

```bash
gmx_mpi trjconv -s cmd.tpr -f cmd.xtc -n gleap.ndx -o cmd.skip10.xtc -skip 10
```



## 7. Summary of Output Directories

| Directory           | Contents                                     |
| ------------------- | -------------------------------------------- |
| `1_modeling/`       | Ligand parameters and force field files      |
| `2_md-preparation/` | System topology, coordinates, index files    |
| `3_min-1/`          | Energy minimized structure                   |
| `4_nvt-1/`          | First NVT equilibration                      |
| `4_nvt-2/`          | Second NVT equilibration                     |
| `5_npt-{1..10}/`    | NPT equilibration with decreasing restraints |
| `6_cmd-{1..5}/`     | Production MD trajectories (5 replicas)      |



## Key Principles

1. **Gradual Restraint Release**: Start with strong position restraints (1000 kJ/mol/nm²) and gradually decrease to 100 kJ/mol/nm² for stable equilibration.

2. **Checkpointing**: Use `.cpt` files to restart interrupted simulations without data loss.

3. **PBC Handling**: Always process PBC artifacts before analysis (whole → cluster → mol).

4. **Multiple Replicas**: Run at least 5 independent simulations for statistical sampling.

5. **Binding Free Energy**: The processed trajectories (COM distances, PCA coordinates) feed into:
   - K-means clustering (20-100 clusters)
   - MSM construction with varying lag times (20-100 frames)
   - PMF calculation using stationary distributions

## Software Requirements

- **AmberTools**: antechamber, tleap, pdb4amber
- **Gaussian**: RESP charge calculation
- **GROMACS**: Molecular dynamics engine
- **Python**: Trajectory analysis and binding free energy calculation
- **deeptime**: Markov state model construction
- **scikit-learn**: K-means clustering



## Citation

When using this workflow, please cite:

- Chen et al., JCIM, 2022 (for the methodology)
- Amber, GROMACS, and deeptime packages as appropriate

