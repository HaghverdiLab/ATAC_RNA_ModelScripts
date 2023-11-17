#!/bin/bash

#$ -N AutoEncHSR
#$ -l h_rt=02:00:00
#$ -cwd  
#$ -o AutoEncModelHSROut.stdlog
#$ -e AutoEncModelHSRErr.stdlog
#$ -l h_vmem=20G
#$ -l m_core=2

source ~/.bashrc
mamba activate babel

MODEL_PATH="/home/sannald/bin/babel_modified"

python3 ${MODEL_PATH}/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/BABEL/HSR_rep7.h5 /fast/AG_Haghverdi/Siddharth_Annaldasula/data/BABEL/HSR_rep8.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_human/hsr
