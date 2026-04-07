# Extract and analyze the top 1000 scoring ligands
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, DataStructs
import matplotlib.pyplot as plt
# import seaborn as sns
import glob
import sys

# in bash
#(head -n 1 out_zinc20_00/affinity.csv && tail -n +2 -q out_zinc20_0{0..3}/affinity.csv | sort -t, -k3,3nr | head -n 1000) > top_1000_affinities.csv
# protein_feats_path_LT=('/chem2/scratch/adlouet/source_code_seb_repo/torch_af_ligand/scripts/inference_ads/protein_feats')

torch_dir = sys.argv[1]
number_ligands_print=sys.argv[2]

protein_paths=glob.glob(torch_dir+'/*/top_1000_affinities.csv')

#protein_paths=glob.glob('/chem2/scratch/adlouet/source_code_seb_repo/torch_af_ligand/scripts/inference_ads/protein_feats/*/top_1000_affinities.csv')
#path='/chem2/scratch/adlouet/source_code_seb_repo/torch_af_ligand/scripts/inference_ads/protein_feats'

#affinity_scores='/chem2/scratch/adlouet/source_code_seb_repo/torch_af_ligand/scripts/inference_ads/protein_feats/pacap_38/top_1000_affinities.csv'

for affinity_scores in protein_paths:
    label = affinity_scores.split('/')[-2]
    individual_pro_path=affinity_scores.split('/top_1000_affinities.csv')[0]
    aff_pre_drop = pd.read_csv(affinity_scores)

    aff_pre_drop['mol_id'] = aff_pre_drop['mol_id'].astype(str).str.upper().str.strip()
    aff = aff_pre_drop.drop_duplicates(subset='mol_id', keep='first')

    zinc = pd.read_csv("/chem2/scratch/adlouet/source_code_seb_repo/torch_af_ligand/scripts/inference_ads/zinc20_combined.csv", sep=' ')
    zinc.columns = [c.strip().lower() for c in zinc.columns]
    zinc['zinc_id'] = zinc['zinc_id'].astype(str).str.upper().str.strip()
    zinc = zinc.drop_duplicates(subset='zinc_id', keep='first')

    merged = aff.merge(zinc, left_on='mol_id', right_on='zinc_id', how='left')
    merged = merged.drop_duplicates(subset='mol_id', keep='first')

    top20 = merged.head(number_ligands_print).dropna(subset=['smiles']).copy()

    top20.to_csv(individual_pro_path+f'/top_{number_ligands_print}_mols.txt')

    mols = [Chem.MolFromSmiles(smi) for smi in top20['smiles']]
    img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(250, 250),
                               legends=[f"{row.mol_id}\n{row.affinity:.2f}" for _, row in top20.iterrows()])

    fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048) for m in mols]
    sim_matrix = [[DataStructs.TanimotoSimilarity(fps[i], fps[j]) for j in range(len(fps))] for i in range(len(fps))]

    sim_df = pd.DataFrame(sim_matrix, index=top20['mol_id'], columns=top20['mol_id'])

    # plt.figure(figsize=(12, 10))
    # sns.heatmap(sim_df, cmap="viridis", annot=True, fmt=".2f", linewidths=0.5, square=True,
    #             cbar_kws={'label': 'Tanimoto Similarity'})
    # plt.xticks(rotation=90)
    # plt.yticks(rotation=0)
    # plt.title(f"Tanimoto Similarity - {label}", fontsize=16, pad=20)
    # plt.tight_layout()
    # plt.show()

    