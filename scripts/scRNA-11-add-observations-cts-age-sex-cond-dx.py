"""
Part 11 of scRNA-seq data analysis. 

Add cell types identified in steps 6-9, age and sex to the adata object.
"""

import scanpy as sc
import numpy as np
import os


# Folder in which result is written
dir_results = './results'


def add_celltype_with_neutrophils(adata):
    """
    Add celltype_with_neutrophils, where neutrophils and MkP2 are reannotated, based on results from parts 6-8.
    """

    observations = adata.obs[['celltype']]
    new_obs = 'celltype_with_neutrophils'
    observations[new_obs] = observations['celltype']

    indices_neutrophils = adata[adata.obs['leiden_myelo-mono_0.08'] == '3'].obs.index
    indices_mkp_2 = adata[adata.obs['leiden_MkP_0.08'] == '1'].obs.index

    observations[new_obs] = observations[new_obs].astype('str')
    observations.loc[indices_neutrophils,new_obs]='Neutrophils'
    observations.loc[indices_mkp_2,new_obs]='Megakaryocyte progenitor 2 (MkP2)'
    observations[new_obs] = observations[new_obs].astype('category')

    adata.obs[new_obs] = observations[new_obs]

    return adata


def add_celltype_short_with_neutrophils(adata):
    """
    Add celltype_short_with_neutrophils, i.e. shorter forms of cell type names.
    """

    dict_Triana_short = {
        'HSCs & MPPs':'HSCs & MPPs',
        'Erythro-myeloid progenitors (EMP)':'EMP',
        'Early erythroid progenitor':'Early eryth prog',
        'Late erythroid progenitor':'Late eryth prog',
        'Aberrant erythroid cells':'Aberrant eryth',
        'Megakaryocyte progenitor (MkP)':'MkP',
        'Megakaryocyte progenitor 2 (MkP2)':'MkP2',
        'Eosinophil/Basophil progenitor (EoBaso)':'EoBaso',
        'Lympho-myeloid prog':'LMP',
        'Early promyelocytes':'Early promyelo',
        'Late promyelocytes':'Late promyelo',
        'Myelocytes':'Myelo',
        'Classical Monocytes':'Class mono',
        'Non-classical monocytes':'Non-class mono',
        'Plasmacytoid dendritic cell progenitors':'pDC prog',
        'Plasmacytoid dendritic cells (pDC)':'pDC',
        'Conventional dendritic cells 1 (cDC1)':'cDC1',
        'Conventional dendritic cells 2 (cDC2)':'cDC2',
        'Pre-pro B cells':'Pre-pro B',
        'Cycling pro-B and pre-B cells':'Cycl pro-B & pre-B',
        'Non-cycling pro-B and pre-B cells':'Non-cycl pro-B & pre-B',
        'Small pre-B cells (light chain re-arrangement)':'Small pre-B',
        'Immature B cells':'Immature B',
        'Mature naïve B cells':'Mature naï B',
        'Non-switched memory B cells':'Non-switch mem B',
        'Class-switched memory B cells':'Class-switch mem B',
        'CD11c+ memory B cells':'CD11c+ mem B',
        'Plasma cells':'Plasma cells',
        'CD4+ naïve T cells':'CD4+ naï T',
        'CD4+ memory T cells':'CD4+ mem T',
        'CD4+ cytotoxic T cells':'CD4+ cytotoxic T',
        'CD4+ CD69+ T cells':'CD4+ CD69+ T',
        'CD8+ naïve T cells':'CD8+ naï T',
        'CD8+ tissue-resident memory T cell':'CD8+ tiss-res mem T',
        'CD8+ central memory T cells ':'CD8+ cent mem T',
        'CD8+ effector memory T cells':'CD8+ eff mem T',
        'gd-T cells':'gd-T',
        'Natural killer T cells':'NK T',
        'CD56dim CD16+ natural killer cells':'CD56dim CD16+ NK',
        'CD56bright CD16- natural killer cells':'CD56bright CD16- NK',
        'Mesenchymal cells 1':'Mesenchymal cells 1',
        'Mesenchymal cells 2':'Mesenchymal cells 2',
        'Natural killer cell progenitor':'NK prog',
        'Immature-like blasts':'Immature blasts',
        'Dendritic-cell like blast':'DC blast',
        'Monocyte-like blast':'Monocyte blast',
        'Neutrophils':'Neutrophils'
    }

    adata.obs['celltype_short'] = adata.obs['celltype'].map(dict_Triana_short)
    adata.obs['celltype_short_with_neutrophils'] = adata.obs['celltype_with_neutrophils'].map(dict_Triana_short)

    return adata


def add_celltype_coarse_with_neutrophils(adata):
    """
    Add celltype_coarse_with_neutrophils, i.e. more coarse cell type categories.
    """

    dict_Triana_coarse = {
        'HSCs & MPPs':'HSCs & MPPs',
        'Erythro-myeloid progenitors (EMP)':'EMP',
        'Early erythroid progenitor':'Erythroid progenitor',
        'Late erythroid progenitor':'Erythroid progenitor',
        'Aberrant erythroid cells':'Aberrant erythroid cells',
        'Megakaryocyte progenitor (MkP)':'MkP',
        'Megakaryocyte progenitor 2 (MkP2)':'MkP',
        'Eosinophil/Basophil progenitor (EoBaso)':'Granulopoietic cells',
        'Lympho-myeloid prog':'LMP',
        'Early promyelocytes':'Granulopoietic cells',
        'Late promyelocytes':'Granulopoietic cells',
        'Myelocytes':'Granulopoietic cells',
        'Classical Monocytes':'Monocytopoietic cells',
        'Non-classical monocytes':'Monocytopoietic cells',
        'Plasmacytoid dendritic cell progenitors':'DC cells',
        'Plasmacytoid dendritic cells (pDC)':'DC cells',
        'Conventional dendritic cells 1 (cDC1)':'DC cells',
        'Conventional dendritic cells 2 (cDC2)':'DC cells',
        'Pre-pro B cells':'B cells',
        'Cycling pro-B and pre-B cells':'B cells',
        'Non-cycling pro-B and pre-B cells':'B cells',
        'Small pre-B cells (light chain re-arrangement)':'B cells',
        'Immature B cells':'B cells',
        'Mature naïve B cells':'B cells',
        'Non-switched memory B cells':'B cells',
        'Class-switched memory B cells':'B cells',
        'CD11c+ memory B cells':'B cells',
        'Plasma cells':'Plasma cells',
        'CD4+ naïve T cells':'T cells',
        'CD4+ memory T cells':'T cells',
        'CD4+ cytotoxic T cells':'T cells',
        'CD4+ CD69+ T cells':'T cells',
        'CD8+ naïve T cells':'T cells',
        'CD8+ tissue-resident memory T cell':'T cells',
        'CD8+ central memory T cells ':'T cells',
        'CD8+ effector memory T cells':'T cells',
        'gd-T cells':'T cells',
        'Natural killer T cells':'NK cells',
        'CD56dim CD16+ natural killer cells':'NK cells',
        'CD56bright CD16- natural killer cells':'NK cells',
        'Mesenchymal cells 1':'Mesenchymal cells',
        'Mesenchymal cells 2':'Mesenchymal cells',
        'Natural killer cell progenitor':'NK cells',
        'Immature-like blasts':'Immature-like blasts',
        'Dendritic-cell like blast':'DC blast',
        'Monocyte-like blast':'Monocyte blast',
        'Neutrophils':'Neutrophils'
    }

    adata.obs['celltype_coarse'] = adata.obs['celltype'].map(dict_Triana_coarse)
    adata.obs['celltype_coarse_with_neutrophils'] = adata.obs['celltype_with_neutrophils'].map(dict_Triana_coarse)
    
    return adata


def add_age(adata):
    """
    Add age at sampling as observation to anndata object, where 0 denotes <18 and 1 denotes >=18
    """

    observations = adata.obs[['sample', 'individual']]

    observations['age_at_sampling'] = ''

    dict_age = {'E1a':'1', 'E1b': '1', 'E1c':'1',
                'E2': '0', 
                'E3': '0', 
                'E4a':'1', 'E4b': '1', 'E4c':'1', 
                'E5': '1', 
                'E6': '1',  
                'E7': '1', 
                'E8a':'1', 'E8b': '1', 'E8c':'1', 
                'E9': '1', 
                'E10':'1', 
                'SDS1': '1', 
                'SDS2': '1', 
                'SDS3': '1', 
                'SDS4': '0', 
                'SDS5a':'1', 'SDS5b':'1'}

    dict_age_2 = {'HCA_BM1':'1', 'HCA_BM2':'1', 'HCA_BM3':'1', 'HCA_BM4':'1', 
                  'HCA_BM5':'1', 'HCA_BM6':'1', 'HCA_BM7':'1', 'HCA_BM8':'1', 
                  'AML1' : '1', 'AML2': '1', 'AML3': '1', 'AML4': '1',
                  'AML5' : '1', 'AML7': '1', 'AML8': '1', 'AML9': '1',
                  'AML10': '1', 'AML11':'1', 'AML12':'1', 'AML13':'1',
                  'AML14': '1', 'AML15':'1', 'AML16':'1', 'AML17':'1',
                  'AML18': '1', 'AML19':'1', 'AML20':'1', 'AML21':'1',
                  'AML_M6':'1'}

    for key, val in dict_age.items():
        observations.loc[observations['sample']==key, 'age_at_sampling'] = val

    for key, val in dict_age_2.items():
        observations.loc[observations['individual']==key, 'age_at_sampling'] = val

    observations['age_at_sampling'] = observations['age_at_sampling'].astype('category')
    adata.obs['age_at_sampling'] = observations['age_at_sampling']

    return adata


def add_sex(adata):
    """
    Add sex as observation to anndata object
    """

    observations = adata.obs[['sample', 'individual']]

    observations['sex'] = ''

    dict_sex = {'E1':'F', 'E2':'M', 'E3':'F', 'E4':'F', 'E5':'F', 
                'E6':'F', 'E7':'M', 'E8':'F', 'E9':'F', 'E10':'M', 
                'SDS1':'M', 'SDS2':'F', 'SDS3':'M', 'SDS4':'M', 'SDS5':'M', 
                'HCA_BM1':'F', 'HCA_BM2':'M', 'HCA_BM3':'M', 'HCA_BM4':'M', 
                'HCA_BM5':'M', 'HCA_BM6':'F', 'HCA_BM7':'F', 'HCA_BM8':'F',
                'AML1':'M', 'AML2':'F', 'AML3':'F', 'AML4':'M', 'AML5':'M', 
                'AML7':'F', 'AML8':'M', 'AML9':'M', 'AML10':'M', 'AML11':'F', 
                'AML12':'M', 'AML13':'F', 'AML14':'M', 'AML15':'M', 'AML16':'M', 
                'AML17':'M', 'AML18':'M', 'AML19':'F', 'AML20':'M', 'AML21':'F', 
                'AML_M6':'M'}

    for key, val in dict_sex.items():
        observations.loc[observations['individual']==key, 'sex'] = val

    observations['sex'] = observations['sex'].astype('category')
    adata.obs['sex'] = observations['sex']

    return adata


def add_condition_short(adata):
    """
    Add condition_short, i.e. a short form of condition, as observation to anndata object
    """

    observations = adata.obs[['condition']]

    new_obs = 'condition_short'
    observations[new_obs] = observations['condition']

    indices_ED = adata[adata.obs['condition'] == 'ERCC6L2 disease'].obs.index
    indices_HD = adata[adata.obs['condition'] == 'healthy donor'].obs.index

    observations[new_obs] = observations[new_obs].astype('str')
    observations.loc[indices_ED,new_obs]='ED'
    observations.loc[indices_HD,new_obs]='HD'
    observations[new_obs] = observations[new_obs].astype('category')

    adata.obs[new_obs] = observations[new_obs]
    return adata


def add_condition_diagnosis(adata):
    """
    Add condition_diagnosis, i.e. the combination of condition_short and diagnosis, as observation to anndata object
    """

    observations = adata.obs[['condition_short','diagnosis']]

    new_obs = 'condition_diagnosis'

    observations['condition_short'] = observations['condition_short'].astype('str')
    observations['diagnosis'] = observations['diagnosis'].astype('str')

    observations[new_obs] = observations['condition_short'] + ' ' + observations['diagnosis']

    indices_HD = adata[adata.obs['condition_short'] == 'HD'].obs.index
    observations.loc[indices_HD,new_obs]='HD'
    indices_AML = adata[adata.obs['diagnosis'] == 'AML'].obs.index
    observations.loc[indices_AML,new_obs]='AML'
    indices_AML_M6 = adata[adata.obs['diagnosis'] == 'AML M6'].obs.index
    observations.loc[indices_AML_M6,new_obs]='AML M6'

    observations[new_obs] = observations[new_obs].astype('category')
    observations['condition_short'] = observations['condition_short'].astype('category')
    observations['diagnosis'] = observations['diagnosis'].astype('category')

    adata.obs[new_obs] = observations[new_obs]

    return adata


def add_condition_diagnosis_BM_normal_separately(adata):
    """
    Add condition_diagnosis_BM_normal_and_BMF, as condition_diagnosis, but BM normal and BMF diagnosis in one group
    """

    observations = adata.obs[['condition_short','diagnosis']]

    new_obs = 'condition_diagnosis_BM_normal_and_BMF'

    observations[new_obs] = observations['condition_diagnosis']
    observations[new_obs] = observations[new_obs].astype('str')

    indices_BM_normal = adata[adata.obs['diagnosis'] == 'BM normal'].obs.index
    observations.loc[indices_BM_normal,new_obs]='ED BM normal and ED BMF'
    indices_BMF = adata[adata.obs['condition_diagnosis'] == 'ED BMF'].obs.index
    observations.loc[indices_BMF,new_obs]='ED BM normal and ED BMF'

    indices_HD = adata[adata.obs['condition_short'] == 'HD'].obs.index
    observations.loc[indices_HD,new_obs]='HD'
    indices_AML = adata[adata.obs['diagnosis'] == 'AML'].obs.index
    observations.loc[indices_AML,new_obs]='AML'
    indices_AML_M6 = adata[adata.obs['diagnosis'] == 'AML M6'].obs.index
    observations.loc[indices_AML_M6,new_obs]='AML M6'

    observations[new_obs] = observations[new_obs].astype('category')

    adata.obs[new_obs] = observations[new_obs]

    return adata


def add_tp53_genotype_wtmut_noTP53clone_wt_HD_wt(adata):
    """
    Add tp53_genotype_wtmut_noTP53clone_wt_HD_wt, where cells are annotated "ED TP53 WT", "ED TP53 MUT", "SDS TP53 WT", "SDS TP53 MUT", and "HD TP53 WT"
    """

    observations = adata.obs[['condition_short','tp53_genotype_wtmut_noTP53clone_wt']]

    new_obs = 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt'

    observations['condition_short'] = observations['condition_short'].astype('str')
    observations['tp53_genotype_wtmut_noTP53clone_wt'] = observations['tp53_genotype_wtmut_noTP53clone_wt'].astype('str')

    observations[new_obs] = observations['condition_short'] + ' TP53 ' + observations['tp53_genotype_wtmut_noTP53clone_wt']

    indices_no_gt = adata[~adata.obs['tp53_genotype_wtmut_noTP53clone_wt'].isin(['WT','MUT'])].obs.index
    observations.loc[indices_no_gt,new_obs]=None

    indices_HD = adata[adata.obs['condition_short'] == 'HD'].obs.index
    observations.loc[indices_HD,new_obs]='HD TP53 WT'

    observations[new_obs] = observations[new_obs].astype('category')

    observations['condition_short'] = observations['condition_short'].astype('category')
    observations['tp53_genotype_wtmut_noTP53clone_wt'] = observations['tp53_genotype_wtmut_noTP53clone_wt'].astype('category')

    adata.obs[new_obs] = observations[new_obs]

    return adata


def add_tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS(adata):
    """
    Add tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS, where cells are annotated "ED BMF TP53 WT", "ED BMF TP53 MUT", "ED MDS TP53 WT", "ED MDS TP53 MUT", etc. and "HD TP53 WT"
    """

    observations = adata.obs[['condition_diagnosis','tp53_genotype_wtmut_noTP53clone_wt']]

    new_obs = 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS'

    observations['condition_diagnosis'] = observations['condition_diagnosis'].astype('str')
    observations['tp53_genotype_wtmut_noTP53clone_wt'] = observations['tp53_genotype_wtmut_noTP53clone_wt'].astype('str')

    observations[new_obs] = observations['condition_diagnosis'] + ' TP53 ' + observations['tp53_genotype_wtmut_noTP53clone_wt']

    indices_no_gt = adata[~adata.obs['tp53_genotype_wtmut_noTP53clone_wt'].isin(['WT','MUT'])].obs.index
    observations.loc[indices_no_gt,new_obs]=None

    indices_HD = adata[adata.obs['condition'] == 'healthy donor'].obs.index
    observations.loc[indices_HD,new_obs]='HD TP53 WT'

    observations[new_obs] = observations[new_obs].astype('category')

    observations['condition_diagnosis'] = observations['condition_diagnosis'].astype('category')
    observations['tp53_genotype_wtmut_noTP53clone_wt'] = observations['tp53_genotype_wtmut_noTP53clone_wt'].astype('category')

    adata.obs[new_obs] = observations[new_obs]

    return adata


# Main

def main():

    adata_file = 'adata_filtered_seeds_scANVI_UMAP_unclear-cell-types_cell-cycle_TP53-gt.h5ad'
    adata = sc.read_h5ad(os.path.join(dir_results, adata_file))

    adata = add_celltype_with_neutrophils(adata)
    adata = add_celltype_short_with_neutrophils(adata)
    adata = add_celltype_coarse_with_neutrophils(adata)
    adata = add_age(adata)
    adata = add_sex(adata)

    adata = add_condition_short(adata)
    adata = add_condition_diagnosis(adata)
    adata = add_condition_diagnosis_BM_normal_separately(adata)

    # TP53 Observations present in anndata object: tp53_genotype_wtmut, tp53_genotype_wtmut_noTP53clone_wt
    adata = add_tp53_genotype_wtmut_noTP53clone_wt_HD_wt(adata)
    adata = add_tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS(adata)


    adata_file_out = 'adata_filtered_seeds_scANVI_UMAP_unclear-cell-types_cell-cycle_TP53-gt_cts-age-sex-cond-dx.h5ad'
    adata.write(os.path.join(dir_results, adata_file_out))


if __name__ == "__main__":
    main()
