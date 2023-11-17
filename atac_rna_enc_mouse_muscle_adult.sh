#!/bin/bash

#$ -N AutoEncMouseMuscleAdult
#$ -l h_rt=24:00:00
#$ -cwd  
#$ -l h_vmem=20G
#$ -o logs/AutoEncModelMouseMuscleAdultBalanced200IntegratedCelltypesout.stdlog
#$ -e logs/AutoEncModelMouseMuscleAdultBalanced200IntegratedCelltypeserr.stdlog

source ~/.bashrc
mamba activate babel

MODEL_PATH="/home/sannald/bin/babel_modified"

python3 /home/sannald/bin/babel_2/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/MouseMuscle_10X/Mouse_Muscle_Adult_ATAC_RNA_Matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_adult/muscle_adult_balanced200_integratedcelltypes_cv --celltype /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_labels/Adult_integrated_celltypes.csv

#python3 /home/sannald/bin/babel_2/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/MouseMuscle_10X/Mouse_Muscle_Adult_ATAC_RNA_Matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_adult/muscle_adult_nocelltypes