"""
Part 12 of scRNA-seq data analysis. 

Normalize scRNA counts.
"""

import os
import scanpy as sc


# Folder in which result is written
dir_results = './results'

def main():

    adata_file = 'adata_filtered_seeds_scANVI_UMAP_unclear-cell-types_cell-cycle_TP53-gt_cts-age-sex-cond-dx.h5ad'
    adata = sc.read_h5ad(os.path.join(dir_results, adata_file))

    sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=True, max_fraction=0.01) # Normalize

    sc.pp.log1p(adata) # Logarithmize the data matrix with X = log(X+1) with natural logarithm

    adata_file_out = 'adata_filtered_seeds_scANVI_UMAP_unclear-cell-types_cell-cycle_TP53-gt_cts-age-sex-cond-dx_norm.h5ad'
    adata.write(os.path.join(dir_results, adata_file_out))


if __name__ == "__main__":
    main()
