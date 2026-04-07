

import os
import pdb
import subprocess
import sys
from tqdm import tqdm


# python cif_to_pdb.py input_files_path output_files_path


# output_pdb='/ptmp/adlouet/mass_produce_md_simulations/protein_simulations'
# input_cif_directory = '/ptmp/adlouet/mass_produce_md_simulations/boltz_output_starting_structures/boltz_results_abeta/predictions/abeta'

output_pdb=str(sys.argv[2])

input_cif_directory = str(sys.argv[1])

# protein_name=input_cif_directory.split('/')[-1]
protein_name = os.path.basename(os.path.normpath(input_cif_directory))

output_pdb_directory=f'{output_pdb}/{protein_name}/apo_protein'
print(input_cif_directory, flush=True)
os.makedirs(output_pdb_directory, exist_ok=True)

cif_files = [f for f in os.listdir(input_cif_directory) if f.endswith(".cif")]
cif_files.sort() 

for cif_file in tqdm(cif_files, desc="Converting files"):
    cif_path = os.path.join(input_cif_directory, cif_file)
    pdb_file = cif_file.replace(".cif", ".pdb")
    pdb_path = os.path.join(output_pdb_directory, pdb_file)

    obabel_command = ["obabel", cif_path, "-O", pdb_path]
    with open(os.devnull, 'w') as devnull:
        subprocess.run(obabel_command, stdout=devnull, stderr=devnull)

    awk_command = ["awk", '/^ATOM |^TER/']
    with open(os.path.join(output_pdb_directory, pdb_file), "rb") as input_file:
        awk_process = subprocess.run(awk_command, stdin=input_file, stdout=subprocess.PIPE, text=True)
        tmp2_content = awk_process.stdout

    with open(os.path.join(output_pdb_directory, pdb_file), "w") as tmp2_file:
        tmp2_file.write(tmp2_content)

print("Conversion complete.")
