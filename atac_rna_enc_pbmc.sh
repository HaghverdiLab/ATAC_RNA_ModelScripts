#!/bin/bash

#$ -N AutoEncPBMCBalanced
#$ -l h_rt=03:00:00
#$ -cwd 
#$ -o logs/AutoEncModelPBMCBalancedFilterOut.stdlog
#$ -e logs/AutoEncModelPBMCBalancedFilterErr.stdlog
#$ -l h_vmem=30G

source ~/.bashrc
mamba activate babel

MODEL_PATH="/home/sannald/bin/babel_modified"

python3 ${MODEL_PATH}/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/BABEL/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_human/pbmc10k_balanced200_filter_cv --celltype /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_human/pbmc10k_celltypes_new.csv
