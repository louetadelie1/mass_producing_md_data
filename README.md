# Mass Producing MD ligand simulations

Scripts for running many protein–ligand molecular dynamics (MD) simulations at once on a GPU cluster (SLURM, MPCDF Viper-GPU).

The target proteins are intrinsically disordered proteins (IDPs), such as Aβ, tau, α-CGRP, PACAP, FUS, TDP-43, p53 and prion. Each one is simulated with a set of top-ranked small molecules from the ZINC library. The goal is to see where, and how strongly, each ligand interacts with the protein.

## What it does

For each protein and each ligand, the pipeline:

1. **Prepares the protein.** It converts Boltz/ColabFold structure predictions (`.cif`) to PDB. Each predicted model becomes a separate replica.
2. **Parameterises the ligand.** It builds GROMACS topologies from SMILES with OpenFF (Sage 2.0).
3. **Builds the system.** It places the ligand next to the protein, adds water and ions, and runs energy minimisation, NVT and NPT equilibration for every replica.
4. **Runs production MD.** All replicas run together on the GPUs with GROMACS (`mdrun -multidir`).
5. **Analyses the results.** It joins the trajectories and computes the distance from the ligand to every residue in every frame.

Everything after step 3 is submitted to SLURM automatically, with each job starting when the previous one finishes.

## Main files

| File | What it does |
|---|---|
| `generate_complex.sh` | Main entry point. Loops over proteins and ligands and submits the jobs. |
| `cif_to_pdb.py` | Converts predicted structures from `.cif` to `.pdb`. |
| `analyze_affinity.py` | Picks the top-scoring ZINC ligands and plots them. |
| `ligand_parameters.py` | Builds ligand force-field files with OpenFF. |
| `run_all_reps_together.sh` | Builds, solvates, minimises and equilibrates each replica. |
| `run_first.sh` | Runs production MD for all replicas. |
| `post_process.sh` | Cleans up and merges the trajectories. |
| `distance_cal.py` | Calculates protein–ligand distances for each residue. |
| `generate_cvs.py` | Writes PLUMED files for metadynamics runs (optional). |
| `run_parameters/` | GROMACS `.mdp` settings: 298 K, 1 bar, 2 fs timestep. |
| `forcefields/` | a99SB-*disp* (default) and CHARMM22* force fields. |
| `previous_cvs_colabfold/` | An older version of the pipeline, kept for reference. |

## Requirements

- SLURM cluster with GPUs
- GROMACS 2025 (MPI build, `gmx_mpi`)
- PLUMED 2 (only for metadynamics)
- Python with `openff-toolkit`, `openff-interchange`, `rdkit`, `mdtraj`, `MDAnalysis`, `numpy`, `pandas`, and `openbabel`

## Usage

1. Edit the paths at the top of `generate_complex.sh` (simulation folder, scripts folder, force field, ligand list).
2. Place the protein replicas at `<simulation_dir>/<protein>/rep_N/rep_N.pdb`.
3. Run:

```bash
bash generate_complex.sh
```

The results for each ligand are written to `<simulation_dir>/<protein>/<ZINC_ID>/`.

## Author

Adélie Louet
