#!/bin/bash

#$ -N AutoEncMouseKidneyAdult
#$ -l h_rt=02:00:00
#$ -cwd  
#$ -o logs/AutoEncModelMouseKidneyAdultOut.stdlog
#$ -e logs/AutoEncModelMouseKidneyAdultErr.stdlog
#$ -l h_vmem=20G

source ~/.bashrc
mamba activate babel

MODEL_PATH="/home/sannald/bin/babel_modified"

python3 ${MODEL_PATH}/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/Mouse_Kidney_10X/M_Kidney_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/kidney_adult 