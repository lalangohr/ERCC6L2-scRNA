"""
Part 3 of scRNA-seq data analysis. 

Identify seeds as input for scANVI.
"""

import scanpy as sc
import os
import numpy as np
import pandas as pd


# Folder in which result is written
dir_results = './results'

# Genes that are problematic for scanpy or scVI
genes_problematic = ['HDD','TRGC','HLA','HLA-DPBI']
    

def parse_marker_list():
    """
    Parse Triana et al. marker genes and prepare them for the identification of seeds
    """

    # Parse marker genes from Triana et al. from their excel
    df_Triana = pd.read_excel('Triana_2021_41590_2021_1059_MOESM3_ESM.xlsx', 
                            sheet_name='Table S4', 
                            header=1, 
                            skiprows=[2,22,24,47])
    df_Triana = df_Triana.rename(columns={'mRNA expression': 'mRNA_expression'})
    df_Triana['mRNA_expression'] = df_Triana.mRNA_expression.apply(lambda x: x.strip(',').split(', '))

    celltype_gene_list = df_Triana['mRNA_expression']

    # Make unique list, i.e. to include each marker gene only once
    celltype_gene_flatlist = [item for sublist in celltype_gene_list for item in sublist]

    celltype_gene_flatlist = [gene.replace('mixture of cells from different AML patients\nCD34', 'CD34') for gene in celltype_gene_flatlist]
    celltype_gene_flatlist = [gene.replace('DNTT expression  can be observed', 'DNTT') for gene in celltype_gene_flatlist]

    gene_set = set(celltype_gene_flatlist) # set of unique genes
    gene_list = list(gene_set) # convert set back to list

    # Remove genes that are problematic for ScanPy or scVI:
    for gene in genes_problematic:
            gene_list.remove(gene)

    # List of curated marker genes is ready
    gene_subset = gene_list

    dict_Triana_pre = {
        'Cycling pro-B and pre-B cells (see text for further information)':'Cycling pro-B and pre-B cells',
        'Non-cycling pro-B and pre-B cells (see text for further information)':'Non-cycling pro-B and pre-B cells',
        'gd-T cells (2 subsets, see text)':'gd-T cells'
    }
    df_Triana = df_Triana.replace({'Cluster label': dict_Triana_pre})


    return df_Triana, gene_subset


def get_score(normalized_adata, gene_set):
    """
    A scANVI helper function.
    Returns the score per cell given a dictionary of + and - genes

    Parameters
    ----------
    normalized_adata
      anndata dataset that has been log normalized and scaled to mean 0, std 1
    gene_set
      a dictionary with two keys: 'positive' and 'negative'
      each key should contain a list of genes
      for each gene in gene_set['positive'], its expression will be added to the score
      for each gene in gene_set['negative'], its expression will be subtracted from its score

    Returns
    -------
    array of length of n_cells containing the score per cell
    """

    score = np.zeros(normalized_adata.n_obs)
    for gene in gene_set['positive']:
        expression = np.array(normalized_adata[:, gene].X)
        score += expression.flatten()
    for gene in gene_set['negative']:
        expression = np.array(normalized_adata[:, gene].X)
        score -= expression.flatten()

    return score


def get_cell_mask(normalized_adata, gene_set):
    """
    A scANVI helper function.
    Calculates the score per cell for a list of genes, then returns a mask for
    the cells with the highest 500 scores. 

    Parameters
    ----------
    normalized_adata
      anndata dataset that has been log normalized and scaled to mean 0, std 1
    gene_set
      a dictionary with two keys: 'positive' and 'negative'
      each key should contain a list of genes
      for each gene in gene_set['positive'], its expression will be added to the score
      for each gene in gene_set['negative'], its expression will be subtracted from its score

    Returns
    -------
    Mask for the cells with the top 50 scores over the entire dataset
    """
    score = get_score(normalized_adata, gene_set)
    cell_idx = score.argsort()[-500:]
    mask = np.zeros(normalized_adata.n_obs)
    mask[cell_idx] = 1

    return mask.astype(bool)


def get_seed_labels(adata, df_Triana, gene_subset):
    """
    Identify seed labels
    """

    # We assign a score to every cell as a function of its cell type signature.
    # In order to compute these scores, we need to normalize the data. 
    # Because this is not the case of scVI or scANVI, we proceed with a copy of the dataset for this step.

    # Normalize data (in a copy of adata)
    normalized = adata.copy()
    sc.pp.normalize_total(normalized, target_sum = 1e4) 

    # Logarithmize the data matrix#
    sc.pp.log1p(normalized)

    # Get subset of data matrix with only marker genes from Triana et al. (or Azimuth et al.)
    normalized = normalized[:,gene_subset].copy()

    # Scale data to unit variance and zero mean
    sc.pp.scale(normalized)

    # Replace cell types names by numbers as scANVI does not handle long cell types properly
    cts = ['HSCs & MPPs', 'Erythro-myeloid progenitors (EMP)', 'Early erythroid progenitor', 'Late erythroid progenitor',
        'Aberrant erythroid cells', 'Megakaryocyte progenitor (MkP)', 'Eosinophil/Basophil progenitor (EoBaso)', 
        'Lympho-myeloid prog', 'Early promyelocytes', 'Late promyelocytes', 'Myelocytes', 'Classical Monocytes',
        'Non-classical monocytes', 'Plasmacytoid dendritic cell progenitors', 'Plasmacytoid dendritic cells (pDC)',
        'Conventional dendritic cells 1 (cDC1)', 'Conventional dendritic cells 2 (cDC2)',
        'Pre-pro B cells', 'Cycling pro-B and pre-B cells', 'Non-cycling pro-B and pre-B cells', 
        'Small pre-B cells (light chain re-arrangement)', 'Immature B cells', 'Mature naïve B cells',
        'Non-switched memory B cells', 'Class-switched memory B cells', 'CD11c+ memory B cells',
        'Plasma cells',
        'CD4+ naïve T cells', 'CD4+ memory T cells', 'CD4+ cytotoxic T cells', 'CD4+ CD69+ T cells',
        'CD8+ naïve T cells', 'CD8+ tissue-resident memory T cell', 'CD8+ central memory T cells ', 
        'CD8+ effector memory T cells', 'gd-T cells',
        'Natural killer T cells', 'CD56dim CD16+ natural killer cells', 'CD56bright CD16- natural killer cells',
        'Mesenchymal cells 1', 'Mesenchymal cells 2',
        'Natural killer cell progenitor',
        'Immature-like blasts', 'Dendritic-cell like blast', 'Monocyte-like blast']
    dict_Triana = {ct: str(idx + 1) for idx, ct in enumerate(cts)}
    df_Triana = df_Triana.replace({'Cluster label': dict_Triana})

    set_to_unknown = False

    # for each cell type and its gene marker list
    for index, row in df_Triana.iloc[::-1].iterrows():
            
            label = row['Cluster label']
            
            curr_gene_list = row['mRNA_expression']
            
            # remove problematic genes 
            for gene in genes_problematic:
                    if gene in curr_gene_list:
                            curr_gene_list.remove(gene)
            
            # remove notes from Triana et al.
            curr_gene_list = [gene.replace('mixture of cells from different AML patients\nCD34', 'CD34') for gene in curr_gene_list]
            curr_gene_list = [gene.replace('DNTT expression  can be observed', 'DNTT') for gene in curr_gene_list]
            
            geneset = {"positive":curr_gene_list,"negative":[]}
            
            print(f'{label} {geneset}')
            
            # call function with normalized subset of adata and gene marker list (dict object) for current cell type
            # returns 50 cells with top score
            mask = get_cell_mask(normalized, geneset,)
            
            if set_to_unknown==False:
                    seed_labels = np.array(mask.shape[0] * ["Unknown"])
                    set_to_unknown = True
            
            seed_labels[mask] = label

    return seed_labels


def main():

    adata_file = os.path.join(dir_results, 'adata_filtered.h5ad')
    adata = sc.read_h5ad(adata_file)

    # Parse Triana et al. marker list
    df_Triana, gene_subset = parse_marker_list()

    # Identify seed labels
    seed_labels = get_seed_labels(adata, df_Triana, gene_subset)

    # Add seed labels to adata object
    adata.obs["seed_labels"] = seed_labels

    # Check what seed label information we have now
    print(adata.obs.seed_labels.value_counts())

    # Save result
    adata.write(os.path.join(dir_results, 'adata_filtered_seeds.h5ad'))

    print('Finished')


if __name__ == "__main__":
    main()
