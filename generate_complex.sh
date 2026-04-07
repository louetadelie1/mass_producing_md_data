#generate_complex_2.sh
# # Need to have already ran LT and analyze_affinity.py
# # need to be in openff environment
# # run analyze_affinity.py

# once the LT is ran, need to have the following dircory setup
# torch_dir='/chem2/scratch/adlouet/source_code_seb_repo/torch_af_ligand/scripts/inference_ads/protein_feats' #LT ouptut directory

"""
dependenceis needed
conda env mass_pro
- openff
- numpy
-  glob
-  sys
-  os
-  pandas
-  MDAnalysis 
-  math
-  shutil
# install pydantic cachetools xmltodict networkx OpenMM pandas numpy glob2 constraint rdkit <-- openff
- colabfold
"""

# Now we have the top ligands — now working with the protein system!
simulation_dir='/ptmp/adlouet/mass_produce_md_simulations/protein_simulations' # where we want to run our simus
scripts='/ptmp/adlouet/mass_produce_md_simulations/scripts'
# fasta_files='/ptmp/adlouet/mass_produce_md_simulations/fasta_files'
# top_molecules='/ptmp/adlouet/mass_produce_md_simulations/top_molecules/zinc_lib/top_50/store_centroid_znc_50.pckl'
run_para='/ptmp/adlouet/mass_produce_md_simulations/scripts/run_parameters'
ff='/ptmp/adlouet/mass_produce_md_simulations/scripts/forcefields/a99SBdisp.ff'
boltz_output_starting_structures='/ptmp/adlouet/mass_produce_md_simulations/boltz_output_starting_structures'
top_molecules_file='/ptmp/adlouet/mass_produce_md_simulations/top_molecules/zinc_lib/top_10/'

cd "${simulation_dir}"

proteins=("alpha_CGRP" "PACAP_27" "PACAP_38" "tau" "medin" "abeta" "fus" "htt" "p53_c" "p53_n" "prion"  "tdp_43")

# this is to generate the protein ensemble -- 
# cd "$simulation_dir" || exit

# for protein_name in "${proteins[@]}"; do
#     echo "Processing: $protein_name"
    
#     #post-processing boltz output
#     # for x in "${boltz_output_starting_structures}"/boltz*; do   for p in ${x}/predictions/*/; do     python "${scripts}"/cif_to_pdb.py "$p" ../protein_simulations/;   done; done

#     dir_name="${protein_name}"

#     cd "$dir_name" || exit
    
#     for i in {0..15}; do
#         mv apo_protein/*_model_${i}.pdb rep_${i}.pdb
#     done

#     for x in *pdb; do
#         mkdir "${x%.*}" && mv "$x" "${x%.*}"
#     done    

#     cd "$simulation_dir" || exit
# done


cd "$simulation_dir"
proteins_="ligand_parameter_file_top_10"
# # this is to generate the ligand parameter files -- 
# #if ligands are for tailored to each protein, must do this seperately:
# for protein_name in "${proteins_[@]}"; do
#     echo "Processing: $protein_name"


#     dir_name="${protein_name}"
#     mkdir -p "$dir_name"
#     cd "$dir_name"
#     current_dir=$(pwd)
#     python $scripts/ligand_parameters.py "$top_molecules_file" "$current_dir"  "$protein_name" 

#     cd "$simulation_dir" || exit
# done

proteins=("abeta")
proteins_="ligand_parameter_file_top_10"

# this is to generate the starting MD ensemble -- 
for protein_name in "${proteins[@]}"; do
    echo "Processing: $protein_name"
    
    cd "$simulation_dir"/"${protein_name}"
    ln -s "$simulation_dir"/"${proteins_}" ./
    mv rep* apo_protein                           
    
    cp -r "$simulation_dir"/"${proteins_}"/* ./

    for x in ZINC*/; do (
        cd "$x" || exit
        echo "$x" 
        dir=$(pwd)

        # cp -r ../apo_protein/* ./
        cp $scripts/run_all_reps_together.sh ./

    #   cp $scripts/parallel_16_reps_prep.sh ./
        sbatch --export=ALL,run_para="$run_para",ff="$ff",protein_name="$protein_name",scripts="$scripts" run_all_reps_together.sh

); done

cd "$simulation_dir" || exit
done


