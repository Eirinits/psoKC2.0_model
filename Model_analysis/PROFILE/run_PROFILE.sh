#!/bin/sh

ulimit -n 3000

model=psoKC

sim_case=updatedModel

list1="Apoptosis,Proliferation,Differentiation,Immune_cells,Inflammation,Th1,Th17,Neutrophils"

i=1
for listnodes in $list1; do
    python3 Scripts/Simulations/MaBoSS_PROFILE_drugs.py $model -sy Linux -p 8 "results_"$sim_case"_WT.txt" -o $listnodes -s "list"$i"_"$sim_case -m "Results/Profiles/WT_stimuli.csv" -rb "Results/Profiles/psoKC_psoriasis_abs_RNA_norm_noOutputs.csv" -rf 100
done