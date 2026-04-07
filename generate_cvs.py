import os
import MDAnalysis as mda
import math
import shutil
import sys


protein_name = sys.argv[1]
path=sys.argv[2]

print("protein_name:", protein_name)
print("path:", path)

pdb_file = os.path.join(path, "template.pdb")

plumed_dir = os.path.join(path, "plumed_files")
os.makedirs(plumed_dir, exist_ok=True)

output_file_cvs = os.path.join(plumed_dir, "CVs.dat")
output_file_cvs_drug = os.path.join(plumed_dir, "CVs_drug.dat")
output_file_master_restart = os.path.join(plumed_dir, "plumed_restart.dat")
output_file_master = os.path.join(plumed_dir, "plumed.dat")

print(pdb_file)
# --- read pdb and extract atom info ---
atoms = []
with open(pdb_file, "r") as f:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            atomnum = int(line[6:11])
            atomname = line[12:16].strip()
            resname = line[17:20].strip()
            resid = int(line[22:26])
            atoms.append((atomnum, atomname, resname, resid))

u = mda.Universe(pdb_file)
number_residues = len(u.residues)-1

# extracting CA/CB atoms - cv3
ca_atoms = [a[0] for a in atoms if a[1] == "CA" and 1 <= a[3] <= number_residues]
cb_atoms = [a[0] for a in atoms if a[1] == "CB" and 1 <= a[3] <= number_residues]

# extracting hydrophobic residues CB - cv4
hydrophobic_residues = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'CYS', 'PRO', 'TRP']
hydro_atoms_cb = [a[0] for a in atoms if a[1] == "CB" and 1 <= a[3] <= number_residues and a[2] in hydrophobic_residues]

# cv5: acidic residues
acid_atoms_ASP = [a[0] for a in atoms if a[2] == "ASP" and a[1] in ("CG", "OD1", "OD2")]
acid_atoms_GLU = [a[0] for a in atoms if a[2] == "GLU" and a[1] in ("CD", "OE1", "OE2")]
acid_atoms = acid_atoms_ASP + acid_atoms_GLU

# basic residues
basic_atoms_ARG = [a[0] for a in atoms if a[2] == "ARG" and a[1] in ("CZ", "NH1", "NH2", "NE")]
basic_atoms_LYS = [a[0] for a in atoms if a[2] == "LYS" and a[1] == "NZ"]
basic_atoms = basic_atoms_ARG + basic_atoms_LYS

# --- write CVs.dat ---
with open(output_file_cvs, "w") as out:
    out.write("#CV1 alpha helix --> inertial\n")
    out.write(f"ALPHARMSD RESIDUES=1-{number_residues} TYPE=DRMSD LESS_THAN={{RATIONAL R_0=0.08 NN=8 MM=12}} LABEL=helix\n\n")

    # out.write("#CV2 sum parallel and anti-parallel beta-content\n")
    # out.write(f"PARABETARMSD RESIDUES=1-{number_residues} TYPE=DRMSD LESS_THAN={{RATIONAL R_0=0.08 NN=8 MM=12}} LABEL=par\n")
    # out.write(f"ANTIBETARMSD RESIDUES=1-{number_residues} TYPE=DRMSD LESS_THAN={{RATIONAL R_0=0.08 NN=8 MM=12}} LABEL=anti\n")
    # out.write("COMBINE ARG=par.lessthan,anti.lessthan PERIODIC=NO LABEL=beta\n\n")

    out.write("#CV3 radius of gyration\n")
    out.write("ca: GROUP ATOMS=" + ",".join(str(x) for x in ca_atoms) + "\n")
    out.write("rgyr: GYRATION ATOMS=ca NOPBC\n\n")

    out.write("#CV4 contacts between hydrophobic residues (CB–CB within 0.6 nm)\n")
    out.write("COORDINATION GROUPA=" + ",".join(str(x) for x in hydro_atoms_cb))
    out.write(" R_0=0.6 NOPBC NLIST NL_CUTOFF=1.0 NL_STRIDE=20 LABEL=hydro\n\n")

    out.write("#CV5 salt bridges (acidic vs basic)\n")
    out.write("COORDINATION GROUPA=" + ",".join(str(x) for x in acid_atoms))
    out.write(" GROUPB=" + ",".join(str(x) for x in basic_atoms))
    out.write(" R_0=0.6 NOPBC NLIST NL_CUTOFF=1.0 NL_STRIDE=20 LABEL=salt\n\n")

    out.write("#CV6 correlation between consecutive psi torsion angles.\n")
    out.write("DIHCOR ...\n")
    for i in range(1, number_residues):
        out.write(f"  ATOMS{i}=@psi-{i},@psi-{i+1}\n")
    out.write("LABEL=dihcor\nNOPBC\n... DIHCOR\n")

# --- write CVs_drug.dat ---
with open(output_file_cvs_drug, "w") as out:
    out.write("#This is the metadynamics CV file for CVs involving the small molecule\n\n")
    out.write("#BELOW ARE SOME CVS BETWEEN 14 regions and the ligand based on CA atoms\n")
    out.write("mol: GROUP NDX_FILE=../../index.ndx NDX_GROUP=MOL0\n\n")

    region = 1
    for i in range(0, len(ca_atoms), 3):
        group = ca_atoms[i:i+3]
        if len(group) < 3:
            break
        start_res = (i // 3) * 3 + 1
        end_res = start_res + 2
        out.write(f"#region {region}, res {start_res}-{end_res}\n")
        out.write("COORDINATION GROUPA=" + ",".join(str(x) for x in group))
        out.write(" GROUPB=mol R_0=1.0 NOPBC NLIST NL_CUTOFF=2.0 NL_STRIDE=20 LABEL=reg" + str(region) + "\n\n")
        region += 1

# --- write plumed.dat ---
protein_atoms = u.select_atoms("protein")
ligand_atoms = u.select_atoms("not protein and not resname SOL and not resname HOH and not name W*")
entity0_start, entity0_end = protein_atoms.atoms[0].id, protein_atoms.atoms[-1].id
entity1_start, entity1_end = ligand_atoms.atoms[0].id, ligand_atoms.atoms[-1].id

n_regions = 0
if os.path.exists(output_file_cvs_drug):
    with open(output_file_cvs_drug) as f:
        for line in f:
            if line.strip().startswith("#region"):
                n_regions += 1
else:
    pass

reg_labels = [f"reg{i}" for i in range(1, n_regions + 1)]
reg_args = ",".join(reg_labels)
num_cvs = 4 + n_regions
bias_factor = round(10 * math.sqrt(num_cvs))
grid_rfiles = [
    "../../rep_0/prod/GRID.helix_lessthan",
    "../../rep_0/prod/GRID.rgyr",
    "../../rep_0/prod/GRID.hydro",
    "../../rep_0/prod/GRID.salt",
    "../../rep_0/prod/GRID.dihcor"
] + [f"../../rep_0/prod/GRID.{r}" for r in reg_labels]


hill_files = [
    "HILLS_helix",
    "HILLS_rgyr",
    "HILLS_hydro",
    "HILLS_salt",
    "HILLS_dihcor"
] + [f"HILLS_{r}" for r in reg_labels]

grid_rfiles_str = ",".join(grid_rfiles)
hill_files_str = ",".join(hill_files)

plumed_text = f"""#this is the master plumed.dat file for metadynamic metainference simulations
RESTART
# define groups
MOLINFO STRUCTURE=../../template.pdb
WHOLEMOLECULES ENTITY0={entity0_start}-{entity0_end} ENTITY1={entity1_start}-{entity1_end}
FLUSH STRIDE=1000

#include file with definition of metadynamics CVs
INCLUDE FILE=../../plumed_files/CVs.dat

#include file with definition of metadynamics CVs involving the small molecule
INCLUDE FILE=../../plumed_files/CVs_drug.dat

#activate parallel bias metadynamics 
PBMETAD ...
    ARG=helix.lessthan,rgyr,hydro,salt,dihcor,{reg_args}
    HEIGHT=1.2
    BIASFACTOR={bias_factor}
    SIGMA=0.64,0.03,0.69,2.75,1.34,{",".join(["1.0"] * n_regions)}
    PACE=500
    GRID_MIN=-1,-1,-1,-1,-1,{",".join(["0"] * n_regions)}
    GRID_MAX=100,20,200,400,50,{",".join(["550"] * n_regions)}
    GRID_SPACING=0.1,0.01,0.1,0.1,0.01,{",".join(["0.1"] * n_regions)}
    GRID_WSTRIDE=10000
    GRID_RFILES={grid_rfiles_str}
    FILE={hill_files_str}
    WALKERS_MPI
    LABEL=pbmetad
... PBMETAD

#print out relevant info
PRINT ARG=helix.lessthan,rgyr,hydro,salt,dihcor,{reg_args},pbmetad.bias FILE=COLVAR STRIDE=500
ENDPLUMED
"""

output_file_master_restart = os.path.join(plumed_dir, "plumed_restart.dat")
with open(output_file_master_restart, "w") as f:
    f.write(plumed_text)

# --- Write the normal plumed.dat (without GRID_RFILES) ---
output_file_master = os.path.join(plumed_dir, "plumed.dat")
lines = []
for line in open(output_file_master_restart):
    if not line.strip().startswith('GRID_RFILES'):
        lines.append(line)

with open(output_file_master, 'w') as file:
    file.writelines(lines)
