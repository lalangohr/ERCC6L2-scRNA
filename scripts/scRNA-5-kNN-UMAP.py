"""
Part 5 of scRNA-seq data analysis. 

Calculate k-NN and UMAP.
"""

import scanpy as sc
import os


# Folder in which result is written
dir_results = './results'


def main():

    adata_file = os.path.join(dir_results, 'adata_filtered_seeds_scANVI.h5ad')
    adata = sc.read_h5ad(adata_file)

    # Compute a neighborhood graph of observations
    sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=30)

    # Embed the neighborhood graph using UMAP (McInnes et al. 2018)
    sc.tl.umap(adata)

    # Save result
    adata.write(os.path.join(dir_results, 'adata_filtered_seeds_scANVI_UMAP.h5ad'))

    print('Finished')


if __name__ == "__main__":
    main()
