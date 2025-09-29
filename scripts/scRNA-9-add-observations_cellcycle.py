"""
Part 9 of scRNA-seq data analysis. 

Calculate cell cycle phase and add it to the adata object.
"""

import scanpy as sc
import numpy as np
import os


# Folder in which result is written
dir_results = './results'


def main():
    adata_file_raw = 'adata_filtered_seeds_scANVI_UMAP_unclear-cell-types.h5ad'
    adata = sc.read_h5ad(os.path.join(dir_results, adata_file_raw))

    adata.layers['counts'] = adata.X.copy() # Save the raw counts in the 'counts' layer

    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True) # Normalize data
    sc.pp.log1p(adata) # Logarithmize the data matrix. 
    sc.pp.scale(adata) # Scale data to unit variance and zero mean.

    # Parse cell cycle genes
    cell_cycle_genes = [x.strip() for x in open('/csc/epitkane/home/llangohr/ERCC6L2-GoT/regev_lab_cell_cycle_genes.txt')]

    # Split into 2 lists
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]

    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

    # Add cell cycle phase by performing cell cycle scoring
    # https://github.com/scverse/scanpy/blob/master/scanpy/tools/_score_genes.py#L191-L267
    # Default phase is S; if G2M is higher than S, it's G2M; if all scores are negative, it's G1.
    sc.tl.score_genes_cell_cycle(
        adata, 
        s_genes=s_genes, 
        g2m_genes=g2m_genes
    )

    adata.X = adata.layers['counts'].copy()    # Put counts back into the .X field before we proceed.

    adata_file_out = 'adata_filtered_seeds_scANVI_UMAP_unclear-cell-types_cell-cycle.h5ad'
    adata.write(os.path.join(dir_results, adata_file_out))


if __name__ == "__main__":
    main()