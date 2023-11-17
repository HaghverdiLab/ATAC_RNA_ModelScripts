#!/bin/bash

#$ -N AutoEncMouseMuscleE14
#$ -l h_rt=02:00:00
#$ -cwd  
#$ -l h_vmem=20G
#$ -o logs/AutoEncModelMouseMuscleE14Balanced200IntegratedCelltypesout.stdlog
#$ -e logs/AutoEncModelMouseMuscleE14Balanced200IntegratedCelltypeserr.stdlog

source ~/.bashrc
mamba activate babel

MODEL_PATH="/home/sannald/bin/babel_modified"

python3 ${BABEL_PATH}/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/MouseMuscle_10X/Mouse_Muscle_E14_ATAC_RNA_Matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_e14/muscle_e14_balanced200_integrated_celltypes_cv --celltype /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_labels/E14_integrated_celltypes.csv

#python3 /home/sannald/bin/babel_2/bin/train_model.py --data /fast/AG_Haghverdi/Siddharth_Annaldasula/data/MouseMuscle_10X/Mouse_Muscle_E14_ATAC_RNA_Matrix.h5 --outdir /fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_e14/muscle_e14_integratedcelltypes
