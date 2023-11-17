import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scanorama
import h5py    
from scipy.io import mmread
from scipy import sparse

# ATAC Regions
interval = pd.read_csv('/fast/AG_Haghverdi/Siddharth_Annaldasula/data/Human_DLPFC_10X/GSE207334_Multiome_atac_peaks.txt.gz', names = ["interval"],compression='gzip',sep='\t', quotechar='"')
interval["interval"] = interval["interval"].apply(lambda x: x.replace("-", ":", 1))

# Values
matrix_rna = mmread('/fast/AG_Haghverdi/Siddharth_Annaldasula/data/Human_DLPFC_10X/GSE207334_Multiome_rna_counts.mtx.gz')
matrix_atac = mmread('/fast/AG_Haghverdi/Siddharth_Annaldasula/data/Human_DLPFC_10X/GSE207334_Multiome_atac_counts.mtx.gz')
matrix = sparse.vstack([matrix_rna,matrix_atac])
matrix_csc = matrix.tocsc()

# Genes
genes = pd.read_csv('/fast/AG_Haghverdi/Siddharth_Annaldasula/data/Human_DLPFC_10X/GSE207334_Multiome_rna_genes.txt.gz', names = ["gene"],compression='gzip',sep='\t', quotechar='"')

# Barcodes
barcodes = pd.read_csv('/fast/AG_Haghverdi/Siddharth_Annaldasula/data/Human_DLPFC_10X/GSE207334_Multiome_cell_meta.txt.gz', compression='gzip', sep='\t', quotechar='"')
barcodes["barcode_only"] = barcodes["barcode"].apply(lambda x: x.split("_")[1])

# Formatting
f = h5py.File('/fast/AG_Haghverdi/Siddharth_Annaldasula/data/Human_DLPFC_10X/Human_DLPFC_Adult_ATAC_RNA_Matrix.h5','w')
grp_matrix = f.create_group("matrix")
grp_matrix.create_dataset("barcodes", data=np.array(barcodes["barcode"], dtype = np.dtype('|S18')))
grp_matrix.create_dataset("data", data=matrix_csc.data)
grp_matrix.create_dataset("indices", data=matrix_csc.indices)
grp_matrix.create_dataset("indptr", data=matrix_csc.indptr)
grp_matrix.create_dataset("shape", data=np.array(matrix_csc.shape))

grp_features = grp_matrix.create_group("features")
grp_features.create_dataset("_all_tag_keys", data=np.array(["genome","interval"], dtype = np.dtype('|S8')))
grp_features.create_dataset("feature_type", data=np.array(features, dtype = np.dtype('|S15')))
grp_features.create_dataset("genome", data=np.full(len(features), "hg38", dtype = np.dtype('|S6')))
grp_features.create_dataset("id", data=np.array(["NA"] * matrix_csc.shape[0], dtype = np.dtype('|S25')))
grp_features.create_dataset("interval", data=np.array((["NA"] * matrix_rna.shape[0]) + list(interval["interval"]), dtype = np.dtype('|S26')))
grp_features.create_dataset("name", data=np.array(list(genes["gene"]) + (["NA"] * matrix_atac.shape[0]), dtype = np.dtype('|S25')))
f.close()