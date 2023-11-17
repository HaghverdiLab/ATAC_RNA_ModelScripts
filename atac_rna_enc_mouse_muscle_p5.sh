#!/bin/bash

#$ -N AutoEncMouseMuscleP5
#$ -l h_rt=24:00:00
#$ -cwd  
#$ -l h_vmem=20G
#$ -o logs/AutoEncModelMouseMuscleP5Balanced200IntegratedCelltypesout.stdlog
#$ -e logs/AutoEncModelMouseMuscleP5Balanced200IntegratedCelltypeserr.stdlog

source ~/.bashrc
mamba activate babel

MODEL_PATH="/home/sannald/bin/babel_modified"

python3 ${MODE_PATH}/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/MouseMuscle_10X/Mouse_Muscle_P5_ATAC_RNA_Matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_p5/muscle_p5_balanced200_integratedcelltypes_cv --celltype /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_labels/P5_integrated_celltypes.csv

#python3 /home/sannald/bin/babel_2/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/MouseMuscle_10X/Mouse_Muscle_P5_ATAC_RNA_Matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_p5/muscle_p5_nocelltypes