#!/bin/bash

#$ -N AutoEncMouseMuscleE18
#$ -l h_rt=02:00:00
#$ -cwd  
#$ -l h_vmem=20G
#$ -o logs/AutoEncModelMouseMuscleE18Balanced200IntegratedCelltypesout.stdlog
#$ -e logs/AutoEncModelMouseMuscleE18Balanced200IntegratedCelltypeserr.stdlog

source ~/.bashrc
mamba activate babel

MODEL_PATH="/home/sannald/bin/babel_modified"

python3 ${MODEL_PATH}/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/MouseMuscle_10X/Mouse_Muscle_E18_ATAC_RNA_Matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_e18/muscle_e18_balanced200_integratedcelltypes_cv --celltype /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_labels/E18_integrated_celltypes.csv 

#python3 /home/sannald/bin/babel_2/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/MouseMuscle_10X/Mouse_Muscle_E18_ATAC_RNA_Matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_e18/muscle_e18_nocelltypes
