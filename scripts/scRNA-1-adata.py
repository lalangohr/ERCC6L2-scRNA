"""
Part 1 of scRNA-seq data analysis. 

Parse data, add metadata and concatenate data for all ERCC6L2 disease and SDS patient, AML, AML M6 and healthy control samples.
"""

import scanpy as sc
import os

print(f'Scanpy version: {sc.__version__}')


# Folder in which result is written. Make sure that dir exists.
dir_results = './results'


def main():
    print('Parse data')

    # Parse data: ERCC6L2 disease
    adata_E1a = sc.read_10x_mtx('E1a/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_E1b = sc.read_10x_mtx('E1b/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_E1c = sc.read_10x_mtx('E1c/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_E2  = sc.read_10x_mtx('E2/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_E3  = sc.read_10x_mtx('E3/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_E4a = sc.read_10x_mtx('E4a/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_E4b = sc.read_10x_mtx('E4b/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_E4c = sc.read_10x_mtx('E4c/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_E5  = sc.read_10x_mtx('E5/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_E6  = sc.read_10x_mtx('E6/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_E7  = sc.read_10x_mtx('E7/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_E8a = sc.read_10x_mtx('E8a/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_E8b = sc.read_10x_mtx('E8b/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_E8c = sc.read_10x_mtx('E8c/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_E9  = sc.read_10x_mtx('E9/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_E10 = sc.read_10x_mtx('E10/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    # Parse data: SDS
    adata_SDS1  = sc.read_10x_mtx('SDS1/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_SDS2  = sc.read_10x_mtx('SDS2/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_SDS3  = sc.read_10x_mtx('SDS3/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_SDS4  = sc.read_10x_mtx('SDS4/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_SDS5a = sc.read_10x_mtx('SDS5a/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_SDS5b = sc.read_10x_mtx('SDS5b/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    # Parse data: AML patients
    adata_AML1  = sc.read_10x_mtx('AML1/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_AML2  = sc.read_10x_mtx('AML2/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_AML3  = sc.read_10x_mtx('AML3/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_AML4  = sc.read_10x_mtx('AML4/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_AML5  = sc.read_10x_mtx('AML5/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_AML7  = sc.read_10x_mtx('AML7/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_AML8  = sc.read_10x_mtx('AML8/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_AML9  = sc.read_10x_mtx('AML9/outs/filtered_feature_bc_matrix/',  var_names='gene_symbols', cache=False)
    adata_AML10 = sc.read_10x_mtx('AML10/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_AML11 = sc.read_10x_mtx('AML11/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_AML12 = sc.read_10x_mtx('AML12/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_AML13 = sc.read_10x_mtx('AML13/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_AML14 = sc.read_10x_mtx('AML14/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_AML15 = sc.read_10x_mtx('AML15/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_AML16 = sc.read_10x_mtx('AML16/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_AML17 = sc.read_10x_mtx('AML17/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_AML18 = sc.read_10x_mtx('AML18/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_AML19 = sc.read_10x_mtx('AML19/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_AML20 = sc.read_10x_mtx('AML20/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',  cache=False)
    adata_AML21 = sc.read_10x_mtx('AML21/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    # Pase AML M6 data
    adata_AML_M6 = sc.read_10x_mtx('AML_M6/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    # Parse HCA data
    adata_HCA_BM1_1 = sc.read_10x_mtx('BM1/cellranger_count_BM1_1/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM1_2 = sc.read_10x_mtx('BM1/cellranger_count_BM1_2/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM1_3 = sc.read_10x_mtx('BM1/cellranger_count_BM1_3/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM1_4 = sc.read_10x_mtx('BM1/cellranger_count_BM1_4/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM1_5 = sc.read_10x_mtx('BM1/cellranger_count_BM1_5/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM1_6 = sc.read_10x_mtx('BM1/cellranger_count_BM1_6/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM1_7 = sc.read_10x_mtx('BM1/cellranger_count_BM1_7/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM1_8 = sc.read_10x_mtx('BM1/cellranger_count_BM1_8/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    adata_HCA_BM2_1 = sc.read_10x_mtx('BM2/cellranger_count_BM2_1/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM2_2 = sc.read_10x_mtx('BM2/cellranger_count_BM2_2/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM2_3 = sc.read_10x_mtx('BM2/cellranger_count_BM2_3/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM2_4 = sc.read_10x_mtx('BM2/cellranger_count_BM2_4/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM2_5 = sc.read_10x_mtx('BM2/cellranger_count_BM2_5/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM2_6 = sc.read_10x_mtx('BM2/cellranger_count_BM2_6/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM2_7 = sc.read_10x_mtx('BM2/cellranger_count_BM2_7/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM2_8 = sc.read_10x_mtx('BM2/cellranger_count_BM2_8/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    adata_HCA_BM3_1 = sc.read_10x_mtx('BM3/cellranger_count_BM3_1/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM3_2 = sc.read_10x_mtx('BM3/cellranger_count_BM3_2/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM3_3 = sc.read_10x_mtx('BM3/cellranger_count_BM3_3/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM3_4 = sc.read_10x_mtx('BM3/cellranger_count_BM3_4/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM3_5 = sc.read_10x_mtx('BM3/cellranger_count_BM3_5/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM3_6 = sc.read_10x_mtx('BM3/cellranger_count_BM3_6/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM3_7 = sc.read_10x_mtx('BM3/cellranger_count_BM3_7/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM3_8 = sc.read_10x_mtx('BM3/cellranger_count_BM3_8/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    adata_HCA_BM4_1 = sc.read_10x_mtx('BM4/cellranger_count_BM4_1/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM4_2 = sc.read_10x_mtx('BM4/cellranger_count_BM4_2/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM4_3 = sc.read_10x_mtx('BM4/cellranger_count_BM4_3/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM4_4 = sc.read_10x_mtx('BM4/cellranger_count_BM4_4/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM4_5 = sc.read_10x_mtx('BM4/cellranger_count_BM4_5/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM4_6 = sc.read_10x_mtx('BM4/cellranger_count_BM4_6/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM4_7 = sc.read_10x_mtx('BM4/cellranger_count_BM4_7/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM4_8 = sc.read_10x_mtx('BM4/cellranger_count_BM4_8/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    adata_HCA_BM5_1 = sc.read_10x_mtx('BM5/cellranger_count_BM5_1/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM5_2 = sc.read_10x_mtx('BM5/cellranger_count_BM5_2/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM5_3 = sc.read_10x_mtx('BM5/cellranger_count_BM5_3/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM5_4 = sc.read_10x_mtx('BM5/cellranger_count_BM5_4/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM5_5 = sc.read_10x_mtx('BM5/cellranger_count_BM5_5/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM5_6 = sc.read_10x_mtx('BM5/cellranger_count_BM5_6/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM5_7 = sc.read_10x_mtx('BM5/cellranger_count_BM5_7/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM5_8 = sc.read_10x_mtx('BM5/cellranger_count_BM5_8/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    adata_HCA_BM6_1 = sc.read_10x_mtx('BM6/cellranger_count_BM6_1/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM6_2 = sc.read_10x_mtx('BM6/cellranger_count_BM6_2/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM6_4 = sc.read_10x_mtx('BM6/cellranger_count_BM6_4/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM6_5 = sc.read_10x_mtx('BM6/cellranger_count_BM6_5/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM6_6 = sc.read_10x_mtx('BM6/cellranger_count_BM6_6/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM6_7 = sc.read_10x_mtx('BM6/cellranger_count_BM6_7/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM6_8 = sc.read_10x_mtx('BM6/cellranger_count_BM6_8/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    adata_HCA_BM7_1 = sc.read_10x_mtx('BM7/cellranger_count_BM7_1/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM7_2 = sc.read_10x_mtx('BM7/cellranger_count_BM7_2/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM7_3 = sc.read_10x_mtx('BM7/cellranger_count_BM7_3/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM7_4 = sc.read_10x_mtx('BM7/cellranger_count_BM7_4/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM7_5 = sc.read_10x_mtx('BM7/cellranger_count_BM7_5/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM7_6 = sc.read_10x_mtx('BM7/cellranger_count_BM7_6/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM7_7 = sc.read_10x_mtx('BM7/cellranger_count_BM7_7/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM7_8 = sc.read_10x_mtx('BM7/cellranger_count_BM7_8/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    adata_HCA_BM8_1 = sc.read_10x_mtx('BM8/cellranger_count_BM8_1/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM8_2 = sc.read_10x_mtx('BM8/cellranger_count_BM8_2/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM8_3 = sc.read_10x_mtx('BM8/cellranger_count_BM8_3/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM8_4 = sc.read_10x_mtx('BM8/cellranger_count_BM8_4/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM8_5 = sc.read_10x_mtx('BM8/cellranger_count_BM8_5/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM8_6 = sc.read_10x_mtx('BM8/cellranger_count_BM8_6/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM8_7 = sc.read_10x_mtx('BM8/cellranger_count_BM8_7/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)
    adata_HCA_BM8_8 = sc.read_10x_mtx('BM8/cellranger_count_BM8_8/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=False)

    # Make list of adata samples
    adata_list = [adata_E1a, adata_E1b, adata_E1c, adata_E2, adata_E3, adata_E4a, adata_E4b, adata_E4c, adata_E5, adata_E6, 
                  adata_E7, adata_E8a, adata_E8b, adata_E8c, adata_E9, adata_E10, 
                  adata_SDS1, adata_SDS2, adata_SDS3, adata_SDS4, adata_SDS5a, adata_SDS5b, 
                  adata_AML1, adata_AML2, adata_AML3, adata_AML4, adata_AML5, adata_AML7, adata_AML8, adata_AML9, adata_AML10, adata_AML11, 
                  adata_AML12, adata_AML13, adata_AML14, adata_AML15, adata_AML16, adata_AML17, adata_AML18, adata_AML19, adata_AML20, adata_AML21,
                  adata_AML_M6,
                  adata_HCA_BM1_1, adata_HCA_BM1_2, adata_HCA_BM1_3, adata_HCA_BM1_4, adata_HCA_BM1_5, adata_HCA_BM1_6, adata_HCA_BM1_7, adata_HCA_BM1_8, 
                  adata_HCA_BM2_1, adata_HCA_BM2_2, adata_HCA_BM2_3, adata_HCA_BM2_4, adata_HCA_BM2_5, adata_HCA_BM2_6, adata_HCA_BM2_7, adata_HCA_BM2_8, 
                  adata_HCA_BM3_1, adata_HCA_BM3_2, adata_HCA_BM3_3, adata_HCA_BM3_4, adata_HCA_BM3_5, adata_HCA_BM3_6, adata_HCA_BM3_7, adata_HCA_BM3_8, 
                  adata_HCA_BM4_1, adata_HCA_BM4_2, adata_HCA_BM4_3, adata_HCA_BM4_4, adata_HCA_BM4_5, adata_HCA_BM4_6, adata_HCA_BM4_7, adata_HCA_BM4_8, 
                  adata_HCA_BM5_1, adata_HCA_BM5_2, adata_HCA_BM5_3, adata_HCA_BM5_4, adata_HCA_BM5_5, adata_HCA_BM5_6, adata_HCA_BM5_7, adata_HCA_BM5_8, 
                  adata_HCA_BM6_1, adata_HCA_BM6_2, adata_HCA_BM6_4, adata_HCA_BM6_5, adata_HCA_BM6_6, adata_HCA_BM6_7, adata_HCA_BM6_8, 
                  adata_HCA_BM7_1, adata_HCA_BM7_2, adata_HCA_BM7_3, adata_HCA_BM7_4, adata_HCA_BM7_5, adata_HCA_BM7_6, adata_HCA_BM7_7, adata_HCA_BM7_8, 
                  adata_HCA_BM8_1, adata_HCA_BM8_2, adata_HCA_BM8_3, adata_HCA_BM8_4, adata_HCA_BM8_5, adata_HCA_BM8_6, adata_HCA_BM8_7, adata_HCA_BM8_8]

    # Metadata
    adata_list_samples = ['E1a', 'E1b', 'E1c', 'E2', 'E3', 'E4a', 'E4b', 'E4c', 'E5', 'E6', 'E7', 'E8a', 'E8b', 'E8c', 'E9', 'E10',
                          'SDS1', 'SDS2', 'SDS3', 'SDS4', 'SDS5a', 'SDS5b', 
                          'AML1', 'AML2', 'AML3', 'AML4', 'AML5', 'AML7', 'AML8', 'AML9', 'AML10', 'AML11', 
                          'AML12', 'AML13', 'AML14', 'AML15', 'AML16', 'AML17', 'AML18', 'AML19', 'AML20', 'AML21',
                          'AML_M6',
                          'HCA_BM1_1', 'HCA_BM1_2', 'HCA_BM1_3', 'HCA_BM1_4', 'HCA_BM1_5', 'HCA_BM1_6', 'HCA_BM1_7', 'HCA_BM1_8',
                          'HCA_BM2_1', 'HCA_BM2_2', 'HCA_BM2_3', 'HCA_BM2_4', 'HCA_BM2_5', 'HCA_BM2_6', 'HCA_BM2_7', 'HCA_BM2_8', 
                          'HCA_BM3_1', 'HCA_BM3_2', 'HCA_BM3_3', 'HCA_BM3_4', 'HCA_BM3_5', 'HCA_BM3_6', 'HCA_BM3_7', 'HCA_BM3_8', 
                          'HCA_BM4_1', 'HCA_BM4_2', 'HCA_BM4_3', 'HCA_BM4_4', 'HCA_BM4_5', 'HCA_BM4_6', 'HCA_BM4_7', 'HCA_BM4_8', 
                          'HCA_BM5_1', 'HCA_BM5_2', 'HCA_BM5_3', 'HCA_BM5_4', 'HCA_BM5_5', 'HCA_BM5_6', 'HCA_BM5_7', 'HCA_BM5_8', 
                          'HCA_BM6_1', 'HCA_BM6_2', 'HCA_BM6_4', 'HCA_BM6_5', 'HCA_BM6_6', 'HCA_BM6_7', 'HCA_BM6_8', 
                          'HCA_BM7_1', 'HCA_BM7_2', 'HCA_BM7_3', 'HCA_BM7_4', 'HCA_BM7_5', 'HCA_BM7_6', 'HCA_BM7_7', 'HCA_BM7_8', 
                          'HCA_BM8_1', 'HCA_BM8_2', 'HCA_BM8_3', 'HCA_BM8_4', 'HCA_BM8_5', 'HCA_BM8_6', 'HCA_BM8_7', 'HCA_BM8_8']
    adata_list_samples_reversed = adata_list_samples.copy()
    adata_list_samples_reversed.reverse()

    adata_list_ind = ['E1', 'E1', 'E1', 'E2', 'E3', 'E4', 'E4', 'E4', 'E5', 'E6', 'E7', 'E8', 'E8', 'E8', 'E9', 'E10', 
                      'SDS1', 'SDS2', 'SDS3', 'SDS4', 'SDS5', 'SDS5', 
                      'AML1', 'AML2', 'AML3', 'AML4', 'AML5', 'AML7', 'AML8', 'AML9', 'AML10', 'AML11', 
                      'AML12', 'AML13', 'AML14', 'AML15', 'AML16', 'AML17', 'AML18', 'AML19', 'AML20', 'AML21',
                      'AML_M6',
                      'HCA_BM1', 'HCA_BM1', 'HCA_BM1', 'HCA_BM1', 'HCA_BM1', 'HCA_BM1', 'HCA_BM1', 'HCA_BM1', 
                      'HCA_BM2', 'HCA_BM2', 'HCA_BM2', 'HCA_BM2', 'HCA_BM2', 'HCA_BM2', 'HCA_BM2', 'HCA_BM2', 
                      'HCA_BM3', 'HCA_BM3', 'HCA_BM3', 'HCA_BM3', 'HCA_BM3', 'HCA_BM3', 'HCA_BM3', 'HCA_BM3', 
                      'HCA_BM4', 'HCA_BM4', 'HCA_BM4', 'HCA_BM4', 'HCA_BM4', 'HCA_BM4', 'HCA_BM4', 'HCA_BM4', 
                      'HCA_BM5', 'HCA_BM5', 'HCA_BM5', 'HCA_BM5', 'HCA_BM5', 'HCA_BM5', 'HCA_BM5', 'HCA_BM5', 
                      'HCA_BM6', 'HCA_BM6', 'HCA_BM6', 'HCA_BM6', 'HCA_BM6', 'HCA_BM6', 'HCA_BM6', 
                      'HCA_BM7', 'HCA_BM7', 'HCA_BM7', 'HCA_BM7', 'HCA_BM7', 'HCA_BM7', 'HCA_BM7', 'HCA_BM7', 
                      'HCA_BM8', 'HCA_BM8', 'HCA_BM8', 'HCA_BM8', 'HCA_BM8', 'HCA_BM8', 'HCA_BM8', 'HCA_BM8']
    adata_list_ind_reversed = adata_list_ind.copy()
    adata_list_ind_reversed.reverse()

    adata_list_condition = ['ERCC6L2 disease','ERCC6L2 disease','ERCC6L2 disease','ERCC6L2 disease','ERCC6L2 disease',
                            'ERCC6L2 disease','ERCC6L2 disease','ERCC6L2 disease','ERCC6L2 disease','ERCC6L2 disease',
                            'ERCC6L2 disease','ERCC6L2 disease','ERCC6L2 disease','ERCC6L2 disease','ERCC6L2 disease',
                            'ERCC6L2 disease',
                            'SDS', 'SDS', 'SDS', 'SDS', 'SDS', 'SDS', 
                            'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 
                            'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 
                            'AML',
                            'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor',
                            'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 
                            'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 
                            'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 
                            'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 
                            'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 
                            'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 
                            'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor', 'healthy donor']
    adata_list_condition_reversed = adata_list_condition.copy()
    adata_list_condition_reversed.reverse()

    adata_list_dx = ['BM normal', 'BM normal', 'BM normal',                         # ERCC6L2 disease
                     'BMF', 'BMF', 'BMF', 'BMF', 'BMF', 'BMF', 'BMF', 'BMF', 'BMF', # ERCC6L2 disease
                     'MDS', 'MDS', 'MDS', 'MDS/AML',                                # ERCC6L2 disease
                     'BMF', 'BMF', 'BMF', 'BMF', 'MDS', 'MDS', # SDS
                     'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', # AML
                     'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', 'AML', # AML
                     'AML M6', # AML M6
                     'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', # HCA
                     'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', # HCA
                     'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', # HCA
                     'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', # HCA
                     'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', # HCA
                     'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy',            # HCA
                     'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', # HCA
                     'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy', 'healthy'] # HCA
    adata_list_dx_reversed = adata_list_dx.copy()
    adata_list_dx_reversed.reverse()

    adata_list_protocol = ['ACK', 'ACK', 'ACK', 'ACK', 'ACK', 'ACK', 'ACK', 'ACK', 'ACK', 'ACK', # ERCC6L2 disease
                           'ACK', 'ACK', 'ACK', 'ACK', 'ACK', 'ACK',                             # ERCC6L2 disease
                           'ACK', 'ACK', 'ACK', 'ACK', 'ACK', 'ACK', # SDS        
                           'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', # AML
                           'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', 'Ficoll', # AML
                           'Ficoll', # AML M6
                           'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', # HCA
                           'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', # HCA
                           'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', # HCA
                           'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', # HCA
                           'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', # HCA
                           'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA',       # HCA
                           'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', # HCA
                           'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'] # HCA
    adata_list_protocol_reversed = adata_list_protocol.copy()
    adata_list_protocol_reversed.reverse()

    adata_list_batch = ['2021', '2021', '2023', '2021', '2021', '2021', '2021', '2023', '2021', '2021', # ERCC6L2 disease
                        '2023', '2021', '2023', '2023', '2023', '2021',                                 # ERCC6L2 disease
                        '2023', '2023', '2023', '2023', '2023', '2023', # SDS
                        '2021', '2021', '2021', '2021', '2021', '2021', '2021', '2021', '2021', '2021', # AML
                        '2021', '2021', '2021', '2021', '2021', '2021', '2021', '2021', '2021', '2021', # AML
                        '2022', # AML M6 - year of publication
                        '2017', '2017', '2017', '2017', '2017', '2017', '2017', '2017', # HCA - year of publication
                        '2017', '2017', '2017', '2017', '2017', '2017', '2017', '2017', # HCA - year of publication
                        '2017', '2017', '2017', '2017', '2017', '2017', '2017', '2017', # HCA - year of publication
                        '2017', '2017', '2017', '2017', '2017', '2017', '2017', '2017', # HCA - year of publication
                        '2017', '2017', '2017', '2017', '2017', '2017', '2017', '2017', # HCA - year of publication
                        '2017', '2017', '2017', '2017', '2017', '2017', '2017',         # HCA - year of publication
                        '2017', '2017', '2017', '2017', '2017', '2017', '2017', '2017', # HCA - year of publication
                        '2017', '2017', '2017', '2017', '2017', '2017', '2017', '2017'] # HCA - year of publication
    adata_list_batch_reversed = adata_list_batch.copy()
    adata_list_batch_reversed.reverse()


    # Add metadata to adata samples
    for adata_curr in adata_list:

        adata_curr.var_names_make_unique()
        
        curr_sample = adata_list_samples_reversed.pop()
        print(curr_sample)
        curr_ind = adata_list_ind_reversed.pop()
        print(curr_ind)
        curr_exp = adata_list_condition_reversed.pop()
        print(curr_exp)
        curr_dx = adata_list_dx_reversed.pop()
        print(curr_dx)
        curr_protocol = adata_list_protocol_reversed.pop()
        print(curr_protocol)
        curr_batch = adata_list_batch_reversed.pop()
        print(curr_batch)
        
        adata_curr.obs['sample'] = curr_sample
        adata_curr.obs['individual'] = curr_ind
        adata_curr.obs['condition'] = curr_exp
        adata_curr.obs['diagnosis'] = curr_dx
        adata_curr.obs['protocol'] = curr_protocol
        adata_curr.obs['year'] = curr_batch # Observation called 'batch' will be overwritten by scVI/scANVI, so add batch as 'year' here
        

    # Concatenate data
    adata = adata_E1a.concatenate(
            adata_E1b, adata_E1c, adata_E2, adata_E3, adata_E4a, adata_E4b, adata_E4c, adata_E5, adata_E6, 
            adata_E7, adata_E8a, adata_E8b, adata_E8c ,adata_E9, adata_E10,
            adata_SDS1, adata_SDS2, adata_SDS3, adata_SDS4, adata_SDS5a, adata_SDS5b,
            adata_AML1, adata_AML2, adata_AML3, adata_AML4, adata_AML5, adata_AML7, adata_AML8, adata_AML9, adata_AML10, adata_AML11, 
            adata_AML12, adata_AML13, adata_AML14, adata_AML15, adata_AML16, adata_AML17, adata_AML18, adata_AML19, adata_AML20, adata_AML21,
            adata_AML_M6,
            adata_HCA_BM1_1, adata_HCA_BM1_2, adata_HCA_BM1_3, adata_HCA_BM1_4, adata_HCA_BM1_5, adata_HCA_BM1_6, adata_HCA_BM1_7, adata_HCA_BM1_8,
            adata_HCA_BM2_1, adata_HCA_BM2_2, adata_HCA_BM2_3, adata_HCA_BM2_4, adata_HCA_BM2_5, adata_HCA_BM2_6, adata_HCA_BM2_7, adata_HCA_BM2_8, 
            adata_HCA_BM3_1, adata_HCA_BM3_2, adata_HCA_BM3_3, adata_HCA_BM3_4, adata_HCA_BM3_5, adata_HCA_BM3_6, adata_HCA_BM3_7, adata_HCA_BM3_8, 
            adata_HCA_BM4_1, adata_HCA_BM4_2, adata_HCA_BM4_3, adata_HCA_BM4_4, adata_HCA_BM4_5, adata_HCA_BM4_6, adata_HCA_BM4_7, adata_HCA_BM4_8, 
            adata_HCA_BM5_1, adata_HCA_BM5_2, adata_HCA_BM5_3, adata_HCA_BM5_4, adata_HCA_BM5_5, adata_HCA_BM5_6, adata_HCA_BM5_7, adata_HCA_BM5_8, 
            adata_HCA_BM6_1, adata_HCA_BM6_2, adata_HCA_BM6_4, adata_HCA_BM6_5, adata_HCA_BM6_6, adata_HCA_BM6_7, adata_HCA_BM6_8, 
            adata_HCA_BM7_1, adata_HCA_BM7_2, adata_HCA_BM7_3, adata_HCA_BM7_4, adata_HCA_BM7_5, adata_HCA_BM7_6, adata_HCA_BM7_7, adata_HCA_BM7_8, 
            adata_HCA_BM8_1, adata_HCA_BM8_2, adata_HCA_BM8_3, adata_HCA_BM8_4, adata_HCA_BM8_5, adata_HCA_BM8_6, adata_HCA_BM8_7, adata_HCA_BM8_8,
            batch_categories=adata_list_samples
        )

    print('adata:')
    print(adata)

    print('adata.obs:')
    print(adata.obs.sample(n=20))

    print('adata.var_names:')
    print(adata.var_names[:5])


    # Save result
    print('Save adata')

    adata.write(os.path.join(dir_results, 'adata.h5ad'))

    print('Finished')



if __name__ == "__main__":
    main()
