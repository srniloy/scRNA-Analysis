

library(Seurat)
library(SeuratDisk)

# rds data

rds_obj <- readRDS('ependymal_cells.rds')


# hdf5 data

hdf5_obj <- Read10X_h5(filename = '20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5',
           use.names = TRUE,
           unique.features = TRUE
           )

seurat_hdf5 <- CreateSeuratObject(counts = hdf5_obj)

# mtx data

mtx_obj <- ReadMtx(mtx = 'raw_feature_bc_matrix/matrix.mtx.gz',
                   cells = 'raw_feature_bc_matrix/barcodes.tsv.gz',
                   features = 'raw_feature_bc_matrix/features.tsv.gz')


seurat_mtx <- CreateSeuratObject(counts = mtx_obj)




