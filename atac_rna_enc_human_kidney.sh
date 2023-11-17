#!/bin/bash

#$ -N AutoEncHumanKidney
#$ -l h_rt=12:00:00
#$ -cwd 
#$ -o logs/AutoEncModelHumanKidneyOut.stdlog
#$ -e logs/AutoEncModelHumanKidneyErr.stdlog
#$ -l h_vmem=30G

source ~/.bashrc
mamba activate babel

MODEL_PATH="/home/sannald/bin/babel_modified"

python3 ${MODEL_PATH}/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/Human_Kidney_10X/H_Kidney_Cancer_Chromium_Nuc_Iso_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_human/human_kidney
