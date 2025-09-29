"""
Part 7 of scRNA-seq data analysis. 

Differential expression (DE) analysis for cell type annotation of unclear cell types, i.e.
DE of cells in MkP and myelo-monocytes cluster compared to cells in other cell types.
"""

import scanpy as sc
import numpy as np
import os


# Folder in which result is written
dir_results = './results'

def main():

    adata_file = 'adata_filtered_seeds_scANVI_UMAP_unclear-cell-types.h5ad'
    adata = sc.read_h5ad(os.path.join(dir_results, adata_file))

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Logarithmize the data matrix. 
    # Computes X=log(X+1), where denotes the natural logarithm unless a different base is given.
    sc.pp.log1p(adata)

    cts_to_cluster = {'MkP':['Megakaryocyte progenitor (MkP)'],
                      'myelo-mono':['Classical Monocytes', 'Early promyelocytes', 'Late promyelocytes', 'Myelocytes']}

    for key, value in cts_to_cluster.items():

        if key == 'MkP':
            clusters=['0','1']
        elif key == 'myelo-mono':
            clusters=['3']

        for cluster in clusters:

            # Get indices of cells from unclear cell type
            indices = adata[adata.obs['leiden_'+key+'_0.08'] ==  cluster].obs.index

            # Split cluster of specfied cell type(s) into own "Unclear" cell type 
            observations = adata.obs[['celltype']]
            ct_key = 'celltype_unclear_'+key
            observations[ct_key] = observations['celltype']
            observations[ct_key] = observations[ct_key].astype('str')
            observations.loc[indices, ct_key]='Unclear'
            observations[ct_key] = observations[ct_key].astype('category')
            adata.obs[ct_key] = observations[ct_key]

            # Remove cell types which Unclear cells were annotated to so that they don't affect Wilcoxon test
            adata_subset = adata[~(adata.obs[ct_key].isin(value))].copy()

            # Rank genes for characterizing groups. Expects logarithmized data
            sc.tl.rank_genes_groups(adata_subset, 
                                    groupby=ct_key,
                                    groups=['Unclear'], 
                                    reference='rest',
                                    method='wilcoxon', 
                                    corr_method='benjamini-hochberg'
            )
            df_markers_2 = sc.get.rank_genes_groups_df(adata_subset, 
                                                       group='Unclear',
                                                       key='rank_genes_groups',
                                                       pval_cutoff=None, 
                                                       log2fc_min=None, 
                                                       log2fc_max=None
            )
            df_markers_2.to_csv(os.path.join(dir_results, 'wilcoxon_'+ct_key+'_cluster'+cluster+'_using_subset.csv'), index=False, sep=',')


if __name__ == "__main__":
    main()
