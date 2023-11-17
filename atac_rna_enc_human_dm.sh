#!/bin/bash

#$ -N AutoEncDM
#$ -l h_rt=04:00:00
#$ -cwd  
#$ -o logs/AutoEncModelDMOut.stdlog
#$ -e logs/AutoEncModelDMErr.stdlog
#$ -l h_vmem=25G
#$ -l m_core=24

source ~/.bashrc
mamba activate babel

MODEL_PATH="/home/sannald/bin/babel_modified"

python3 ${MODEL_PATH}/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/BABEL/DM_rep4.h5 /fast/AG_Haghverdi/Siddharth_Annaldasula/data/BABEL/DM_rep8.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_human/dm