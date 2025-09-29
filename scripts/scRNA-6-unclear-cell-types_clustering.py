"""
Part 6 of scRNA-seq data analysis. 

Cluster MkP and myelo-monocytes.
"""

import scanpy as sc
import numpy as np
import os


# Folder in which result is written
dir_results = './results'

def main():

    adata_file = 'adata_filtered_seeds_scANVI_UMAP.h5ad'
    adata = sc.read_h5ad(os.path.join(dir_results, adata_file))


    cts_to_cluster = {'MkP':['Megakaryocyte progenitor (MkP)'],
                    'myelo-mono':['Classical Monocytes', 'Early promyelocytes', 'Late promyelocytes', 'Myelocytes']}

    for key, value in cts_to_cluster.items():

        # Cluster cells of specified cell types
        adata_to_cluster = adata[adata.obs['celltype'].isin(value)]

        # Normalize
        sc.pp.normalize_total(adata_to_cluster, target_sum=1e4)

        # Logarithmize the data matrix. 
        # Computes X=log(X+1), where denotes the natural logarithm unless a different base is given.
        sc.pp.log1p(adata_to_cluster)

        # Scale data to unit variance and zero mean.
        sc.pp.scale(adata_to_cluster)

        resolutions = [0.06, 0.08, 0.1, 0.2, 0.5, 1.0]  # default resolution in 1.0
        for resolution in resolutions:
            key_to_add = 'leiden_'+key+'_'+str(resolution)
            sc.tl.leiden(adata_to_cluster, resolution = resolution, key_added = key_to_add)
            adata.obs[key_to_add] = adata_to_cluster.obs[key_to_add]

    # Save file
    adata_file_out = 'adata_filtered_seeds_scANVI_UMAP_unclear-cell-types.h5ad'
    adata.write(os.path.join(dir_results, adata_file_out))


if __name__ == "__main__":
    main()
