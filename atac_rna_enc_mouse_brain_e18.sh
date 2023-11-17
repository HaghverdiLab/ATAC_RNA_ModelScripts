#!/bin/bash

#$ -N AutoEncMouseBrainE18
#$ -l h_rt=03:00:00
#$ -cwd  
#$ -o logs/AutoEncModelMouseE18BrainOut.stdlog
#$ -e logs/AutoEncModelMouseE18BrainErr.stdlog
#$ -l h_vmem=20G

source ~/.bashrc
mamba activate babel

MODEL_PATH="/home/sannald/bin/babel_modified"

python3 ${MODEL_PATH}/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/MouseBrain_10X/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/brain_e18