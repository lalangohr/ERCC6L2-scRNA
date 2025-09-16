"""
Part 2 of scRNA-seq data analysis. 

Filter low quality cells, filter dying cells, i.e. cells with high amount of mtRNA transcript.
"""

import scanpy as sc
import os


# Folder in which result is written
dir_results = './results'


def main():

    adata_file = os.path.join(dir_results, 'adata.h5ad')
    adata = sc.read_h5ad(adata_file)    
        
    # Filtering low quality cells
    sc.pp.filter_genes(adata, min_cells=3)   # Filter genes: only keep genes with at least min_counts=3
    sc.pp.filter_cells(adata, min_genes=200) # Filter cells: only keep cells with at least min_genes=200
         
    # Annotate mt, ribo and hem genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")            # Mitochondrial genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL")) # Ribosomal genes
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))      # Hemoglobin genes
        
    # Calculate quality control metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
        
    # Filtering cells with high amount of mtRNA transcripts 
    adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > 20
    adata.obs.mt_outlier.value_counts()
    adata = adata[(~adata.obs.mt_outlier)].copy()
        
    # Save result
    adata.write(os.path.join(dir_results, 'adata_filtered.h5ad'))
    print('Finished')


if __name__ == "__main__":
    main()
