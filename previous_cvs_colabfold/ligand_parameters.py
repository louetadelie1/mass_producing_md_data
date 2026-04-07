from openff.units import unit
import numpy as np
from openff.interchange import Interchange
from openff.toolkit import ForceField
from openff.toolkit import Molecule, Topology
import glob
import sys
import os
import pandas as pd

# protein_feats_path_LT=('/chem2/scratch/adlouet/source_code_seb_repo/torch_af_ligand/scripts/inference_ads/protein_feats')

# protein_paths=glob.glob(f'/chem2/scratch/adlouet/source_code_seb_repo/torch_af_ligand/scripts/inference_ads/protein_feats/{$base}/top_20_mols.txt')


base = sys.argv[1]
current_dir=sys.argv[2]
protein_name=sys.argv[3]
top_mols=sys.argv[4]

# base = 'p53_n'
# current_dir = '/chem2/scratch/adlouet/mass_produce_md_simulations/protein_simulations/p53_n'

# base = 'p53_n'
# current_dir = '/chem2/scratch/adlouet/mass_produce_md_simulations/protein_simulations/p53_n'

# top_molecules = f"/chem2/scratch/adlouet/source_code_seb_repo/torch_af_ligand/scripts/inference_ads/protein_feats/{protein_name}/top_20_mols.txt"
top_molecules = f"{top_mols}/top_20_mols_{protein_name}.txt"

sage = ForceField("openff-2.0.0.offxml")
cubic_box = unit.Quantity(30 * np.eye(3), unit.angstrom)

df = pd.read_csv(top_molecules)

for i, row in df.iterrows():
    smiles = row['smiles']
    zinc_id = str(row['mol_id'])

    try:
        mol = Molecule.from_smiles(smiles)
        mol.generate_conformers()

        interchange = Interchange.from_smirnoff(
            topology=[mol],
            force_field=sage,
            box=cubic_box
        )
        zinc_dir = os.path.join(current_dir, zinc_id)
        param_dir = os.path.join(zinc_dir, f"{zinc_id}.param")

        os.makedirs(param_dir, exist_ok=True)

        output_prefix = os.path.join(param_dir, "system")
        # output_dir = os.path.join(current_dir, zinc_id)
        # os.makedirs(output_dir, exist_ok=True)

        # output_prefix = os.path.join(output_dir, "system")
        interchange.to_gromacs(output_prefix, monolithic=False)
        print(f'Wrapped up',zinc_id)
    except Exception as e:
        print(f" {smiles} didn't work ({e})")

