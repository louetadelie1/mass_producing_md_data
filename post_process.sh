printf "1 | 13\n1\nq\n" | gmx_mpi make_ndx -f rep_0/prod/prod_final.tpr

for y in rep*/; do
  cd "${y}/prod" || continue
  echo "  Concatenating parts in $x/$y/prod ..."
  gmx_mpi trjcat -f traj_comp.part* -o concatenated.xtc -cat
  echo 20 | gmx_mpi trjconv -f concatenated.xtc  -s *tpr -pbc mol -ur compact -o pbc_mol.xtc -n ../../index.ndx
  cd ../..   # back to ph*/ directory
done

# Step 2: concatenate all replica trajectories into one file
echo "Combining all replicas for $x ..."
gmx_mpi trjcat -f rep*/prod/pbc_mol.xtc -o all_concatenated.xtc -cat
echo 20 | gmx_mpi editconf -f rep_0/prod/prod_final.tpr -o template.pdb -n index.ndx 

PWD_VAR=$(pwd)
python distance_cal.py PWD_VAR="$PWD_VAR"

