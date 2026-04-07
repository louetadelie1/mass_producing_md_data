#!/bin/bash -l
#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J zn_run_prep
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4

module purge
module load gcc/13 cmake/3.30 openmpi/5.0
module load gromacs/2025.1

plumed_prefix=/u/adlouet/src/plumed2
export PLUMED_KERNEL=$plumed_prefix/lib/libplumedKernel.so
export LD_LIBRARY_PATH=$plumed_prefix/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export GMX_FORCE_UPDATE_DEFAULT_GPU=1
export GMX_DISABLE_GPU_TIMING=1



# Number of concurrent jobs to run
MAX_JOBS=4   # You requested 4 tasks per node; so allow 4 replicas in parallel
JOB_COUNT=0

for rep in rep*; do

  (
    echo "=== Processing $rep ==="
    cd "$rep" || exit

    # -------------------------------
    # Your whole replica pipeline
    # Everything stays the same
    # -------------------------------
    mkdir -p generate_complex
    cd generate_complex || exit
    ln -s "$ff" ./

    gmx_mpi pdb2gmx -f ../*.pdb -o "$protein_name.gro" -p "$protein_name.top" -ignh <<EOF
1
1
EOF


    ln -s ../../*param ./

    cp *param/*_MOL0.itp Ligand.itp
    cp *param/*_atomtypes.itp Ligand_atomtypes.itp
    cp *param/system.gro liga.gro

      cp "$protein_name.top" Complex2.top

    sed -e '/#include ".\/a99SBdisp.ff\/forcefield.itp"/a\
\
; Include ligand parameters\
#include "Ligand_atomtypes.itp"\
\
' -e '/; Include water topology/i\
\
; Include ligand topology\
#include "Ligand.itp"\
\
' Complex2.top > Complex.top
    echo "MOL0  1" >> Complex.top

    # ln -s ../../*param ./

    # cp *param/*_MOL0.itp Ligand.itp
    # cp *param/*_atomtypes.itp Ligand_atomtypes.itp
    # cp *param/system.gro liga.gro

    # cp "$protein_name.top" Complex2.top

    # sed -e '/#include ".\/a99SBdisp.ff\/forcefield.itp"/a\
    # ; Include ligand parameters\
    # #include "Ligand_atomtypes.itp"\
    # ' -e '/; Include water topology/i\
    # ; Include ligand topology\
    # #include "Ligand.itp"\
    # ' Complex2.top > Complex.top
    # echo "MOL0  1" >> Complex.top

    gmx_mpi editconf -f "$protein_name.gro" -o tiny_apo_pro_box.gro -c -bt cubic -d 0.1
    gmx_mpi insert-molecules -f tiny_apo_pro_box.gro -ci liga.gro -nmol 1 -radius 0.05 -try 50000 -o Complex.gro

    cd ../
    mkdir -p prep
    cd prep || exit
    ln -s "$ff" ./

    gmx_mpi editconf -f ../generate_complex/Complex.gro -o Complex_box.gro -c -d 1.5 -bt dodecahedron
    gmx_mpi solvate -cp Complex_box.gro -cs tip4p.gro -o Complex_solv.gro -p ../generate_complex/Complex.top
    gmx_mpi grompp -f "$run_para/em_ions.mdp" -c Complex_solv.gro -p ../generate_complex/Complex.top -o Complex_solv.tpr -maxwarn 3

    grep "qtot" ../generate_complex/Complex.top | tail -n 1 | awk '{print $11}'

    gmx_mpi genion -s Complex_solv.tpr -o Complex_ionized.gro -p ../generate_complex/Complex.top -nname CL -pname NA -neutral <<EOF
15
EOF

    cd ../
    mkdir -p minim/em
    cd minim/em || exit

    gmx_mpi grompp -f "$run_para/minim.mdp" -c ../../prep/Complex_ionized.gro -p ../../generate_complex/Complex.top -o Complex_em.tpr -maxwarn 1
    gmx_mpi mdrun -v -s Complex_em.tpr

    cd ../../
    mkdir -p minim/nvt
    cd minim/nvt || exit

    gmx_mpi make_ndx -f ../../generate_complex/Complex.gro <<EOF
1
13
q
EOF

    gmx_mpi genrestr -f ../../generate_complex/liga.gro -n index.ndx -o ../../generate_complex/posre_LIGAND.itp -fc 1000 1000 1000 <<EOF
13
EOF

    sed -i '/#include "ligand.itp"/a \
    ; Ligand position restraints\n#ifdef POSRES\n#include "posre_LIGAND.itp"\n#endif' ../../generate_complex/Complex.top

    gmx_mpi grompp -f "$run_para/nvt.mdp" -c ../em/confout.gro -r ../em/confout.gro -p ../../generate_complex/Complex.top -o Complex_nvt.tpr -maxwarn 1
    gmx_mpi mdrun -v -s Complex_nvt.tpr

    cd ../
    mkdir -p npt
    cd npt || exit

    gmx_mpi grompp -f "$run_para/npt.mdp" -c ../nvt/confout.gro -r ../nvt/confout.gro -p ../../generate_complex/Complex.top -o Complex_npt.tpr -maxwarn 15
    gmx_mpi mdrun -v -s Complex_npt.tpr

    cd ../../
    mkdir -p prod
    cd prod || exit

    gmx_mpi grompp -f "$run_para/run.mdp" -c ../minim/npt/confout.gro -r ../minim/npt/confout.gro -p ../generate_complex/Complex.top -o prod_final.tpr -maxwarn 15

    echo "=== Finished $rep ==="
  ) &   # <<<<<<<<<<<<<<<< BACKGROUND JOB

  ((JOB_COUNT++))

  # Limit concurrency to the number of CPU tasks per node (4)
  if [[ $JOB_COUNT -ge $MAX_JOBS ]]; then
      wait
      JOB_COUNT=0
  fi

done

wait   # Make sure all background jobs finished before continuing


cp *1_seed_001.pdb/generate_complex/Complex.gro template.gro
cp *1_seed_001.pdb/minim/nvt/index.ndx ./    
gmx_mpi editconf -f template.gro -o template.pdb
sed -i 's/OC1/  O/' template.pdb
sed -i 's/OC2/OXT/' template.pdb

base_ZN=$(basename "$PWD")
path_=$simulation_dir/$protein_name/$base_ZN
rm -r plumed_files
python  $scripts/generate_cvs_files.py "$path_" "$protein_name"

i=0
for d in rep_*; do
  mv "$d" "rep_${i}"
  ((i++))
done

rm -f slurm*

for x in rep*; do (cd ${x}/prod && rm -f !(*tpr));done

cp "$run_sh" ./
cp "$run_restart" ./

sed -i "s/protein/${protein_name}/g" run_first.sh
sed -i "s/protein/${protein_name}/g" run_restart.sh

jobid=$(sbatch --nice 100 run_first.sh | awk '{print $4}')
sbatch --dependency=afterok:$jobid run_restart.sh