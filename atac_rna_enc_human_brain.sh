#!/bin/bash

#$ -N AutoEncBrain
#$ -l h_rt=03:00:00
#$ -cwd  
#$ -o logs/AutoEncModelBrainOut.stdlog
#$ -e logs/AutoEncModelBrainErr.stdlog
#$ -l h_vmem=20G

source ~/.bashrc
mamba activate babel

MODEL_PATH="/home/sannald/bin/babel_modified"

python3 ${MODEL_PATH}/bin/train_model.py --data /scratch/AG_Haghverdi/Siddharth_Annaldasula/data/Brain_Human10XMulti/human_brain_3k_filtered_feature_bc_matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_human/brain