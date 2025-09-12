# gene2gene script to calculate alignment of pseudotime trajectories
# pseudotime calculated with celloracle

import warnings
import scanpy as sc
# import palantir
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sb
from pathlib import Path
import igraph as ig
import sys
import anndata as ad


# celloracle
import celloracle as co
from celloracle.applications import Pseudotime_calculator

# g2g
from genes2genes import Main
from genes2genes import MyFunctions 
from genes2genes import TimeSeriesPreprocessor
from genes2genes import PathwayAnalyserV2
from genes2genes import VisualUtils
from genes2genes import ClusterUtils
from leven import levenshtein  
from sklearn.metrics.pairwise import euclidean_distances
from scipy.spatial import distance
from tqdm import tqdm
import scipy
from sklearn.cluster import AgglomerativeClustering
import sklearn 
import matplotlib.patches as mpatches
from pathlib import Path
from optbinning import ContinuousOptimalBinning
from scipy import stats
import shapely
from shapely.geometry import LineString
from matplotlib.lines import Line2D

ROOT_ID = ""
CELL_CLUSTER = "celltype_short_with_neutrophils"
JOINT_CMAP={'HSCs & MPPs':'#029E73', 'EMP':'#DE8F05' , 'Early eryth prog':'#0173B2' , 'Late eryth prog':'#D55E00', 
    'Mature naï B': '#0173B2','Non-switch mem B': '#DE8F05','Cycl pro-B & pre-B': '#029E73',
    'Non-cycl pro-B & pre-B': '#D55E00','Small pre-B': '#CC78BC','LMP': '#CA9161','Immature B': '#FBAFE4'}

## CATEGORIES
HD_CASES = ["HD","HD TP53 WT"]
HD_SELECTION_COL = ['condition_diagnosis', 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt']
ED_CASES = ["ED BMF", "ED MDS","ED MDS/AML","ED MDS/AML-clus1","ED TP53 WT", "ED TP53 MUT",'ED BMF TP53 WT', 'ED MDS TP53 WT', 'ED MDS TP53 MUT']
ED_SELECTION_COL = ['condition_diagnosis', 'condition_diagnosis', 'condition_diagnosis', 'condition_diagnosis',
                    'tp53_genotype_wtmut_noTP53clone_wt_HD_wt', 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt',
                    'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS', 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS', 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS']
SDS_CASES = ["SDS BMF", "SDS MDS", "SDS TP53 WT", "SDS TP53 MUT", 'SDS BMF TP53 WT', 'SDS MDS TP53 WT', 'SDS MDS TP53 MUT']
SDS_SELECTION_COL = ['condition_diagnosis', 'condition_diagnosis', 
                    'tp53_genotype_wtmut_noTP53clone_wt_HD_wt', 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt',
                    'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS', 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS', 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS']

## CELLTYPES
CT_ERY = ['HSCs & MPPs', 'EMP', 'Early eryth prog', 'Late eryth prog']
CT_THROMBO = ['HSCs & MPPs', 'EMP', 'MkP'] 
CT_B = ['HSCs & MPPs','LMP', 'Pre-pro B', 'Cycl pro-B & pre-B', 'Non-cycl pro-B & pre-B', 'Immature B', 'Mature naï B', 'Non-switch mem B', 'Small pre-B']
CT_T_NK = ['HSCs & MPPs','LMP','CD4+ cytotoxic T', 'CD4+ mem T', 'CD4+ naï T', 'CD8+ cent mem T', 'CD8+ naï T', 'CD56bright CD16 NK', 'NK T']
CT_DENDRITIC = ['HSCs & MPPs', 'EMP','pDC prog', 'cDC1']
CT_MYELO = ['HSCs & MPPs', 'EMP','Early promyelo', 'Late promyelo', 'Myelo', 'Class mono', 'Non-class mono']
LINEAGE_DICT = {"Lineage_Eryth": CT_ERY,
                "Lineage_Thrombo": CT_THROMBO,
                "Lineage_B": CT_B,
                "Lineage_T_Nk": CT_T_NK,
                "Lineage_Dendritic": CT_DENDRITIC,
                "Lineage_Myelo": CT_MYELO,
                }

ERY_MARKER = ['CD34', 'GATA2', 'KIT',  'TFRC', 'EPOR', 'CD36', 'GYPA', 'SLC4A1', 'MYC', 'MYB', 'ITGA4','GATA1', 'SBDS'] 
B_MARKER = ['IKZF1', 'SPI1', 'MYC', 'MYB', 'RUNX1', 'CBFB', 'TCF3', 'EBF1', 'PAX5', 'MME', 'CD19', 'BCL11A', 'PTPRC', 'SELL']
GERM_SOM_MARKER = ['ERCC6L2', 'SBDS', 'TP53']
LEAD_MARKER = ['MYC', 'BCL2L1', 'SON', 'SOD2', 'RPS14', 'EZH2', 'TFRC', 'CD36', 'ITGA4', 'GCLM']

## PT COMPARISONS
# HD vs ED 
HD_REF_ED = [
    "HD", "HD", "HD", "ED BMF", "ED BM normal", "ED BM normal",#"HD", "ED BMF",#condition diagnosis
    "HD TP53 WT", "HD TP53 WT", "ED TP53 WT", #tp53_genotype_wtmut_noTP53clone_wt_HD_wt
    "HD TP53 WT", 'ED BMF TP53 WT', #"HD TP53 WT", "HD TP53 WT",'ED MDS TP53 WT',  # 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS'
]
ED_QUERY_HD =[
    "ED BM normal","ED BMF", "ED MDS/AML", "ED MDS/AML", "ED MDS/AML","ED BMF", # "ED BMF""ED MDS/AML-clus1", "ED MDS/AML-clus1", ,#condition diagnosis
    "ED TP53 WT", "ED TP53 MUT", "ED TP53 MUT", #tp53_genotype_wtmut_noTP53clone_wt_HD_wt
    "ED BMF TP53 WT", "ED MDS TP53 WT", #"ED MDS TP53 MUT",'ED MDS TP53 MUT', 'ED MDS TP53 WT' #'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS'
    ]
HD_ED_COL = ['condition_diagnosis' for i in range(6)] 
HD_ED_COL += ['tp53_genotype_wtmut_noTP53clone_wt_HD_wt' for i in range(3)] 
HD_ED_COL += ['tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS' for i in range(2)]

# HD vs SDS
HD_REF_SDS = [
        "HD", "HD", "SDS BMF", #condition diagnosis
        "HD TP53 WT", "HD TP53 WT", "SDS TP53 WT", #tp53_genotype_wtmut_noTP53clone_wt_HD_wt
        "HD TP53 WT", "HD TP53 WT", "HD TP53 WT",'SDS MDS TP53 WT', 'SDS BMF TP53 WT' # 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS'
    ]
SDS_QUERY_HD = [
        "SDS BMF", "SDS MDS", "SDS MDS",#condition diagnosis
        "SDS TP53 WT", "SDS TP53 MUT", "SDS TP53 MUT", #tp53_genotype_wtmut_noTP53clone_wt_HD_wt
        "SDS BMF TP53 WT", "SDS MDS TP53 WT", "SDS MDS TP53 MUT",'SDS MDS TP53 MUT', 'SDS MDS TP53 WT' #'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS'
        ]
HD_SDS_COL = ['condition_diagnosis' for i in range(3)] 
HD_SDS_COL += ['tp53_genotype_wtmut_noTP53clone_wt_HD_wt' for i in range(3)] 
HD_SDS_COL += ['tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS' for i in range(5)]

# ED vs SDS
ED_REF_SDS = [
        "ED BMF", "ED MDS", "ED MDS/AML", #condition diagnosis
        "ED TP53 WT", "ED TP53 MUT", #tp53_genotype_wtmut_noTP53clone_wt_HD_wt
        "ED BMF TP53 WT", #"ED MDS TP53 WT", "ED MDS TP53 MUT" #'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS'
        ]
SDS_QUERY_ED = [
        "SDS BMF", "SDS MDS", "SDS MDS",#condition diagnosis
        "SDS TP53 WT", "SDS TP53 MUT",  #tp53_genotype_wtmut_noTP53clone_wt_HD_wt
        "SDS BMF TP53 WT", #'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS'
        ]
ED_SDS_COL = ['condition_diagnosis' for i in range(3)] 
ED_SDS_COL += ['tp53_genotype_wtmut_noTP53clone_wt_HD_wt' for i in range(2)] 
ED_SDS_COL += ['tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS']


## UTILS 
def remove_suffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[:-len(suffix)]
    return input_string

def merge_files(filename1, filename2):
    df1 = pd.read_csv(filename1, index_col=0) 
    df2 = pd.read_csv(filename2, index_col=0) 
    merged_df = pd.concat([df1, df2])
    print(merged_df)
    
    name = remove_suffix(filename1,'.csv')
    print(name)
    merged_df.to_csv(name+'merged.csv')

    return


## Visualization PT Alignment
# visualisations from g2g that are not accessible yet in the git repository, copied code from gene2gene repo
def plot_alignmentSim_vs_l2fc(df, path_fig):
    plt.subplots(1,2,figsize=(16,6))
    # histogram
    plt.subplot(1,2,1)
    temp = np.asarray([x*100 for x in df['alignment_similarity_percentage']]) 
    p = sb.histplot(temp, label='Alignment Similarity %', bins=10)
    plt.xlabel('Alignment Similarity Percentage')
    plt.xlim([0,100])
    p.set_yticklabels(p.get_yticks(), size = 12)
    p.set_xticklabels(p.get_xticks(), size = 12)
    plt.xlabel('Alignment similarity percentage', fontsize='12')
    plt.ylabel('Frequency', fontsize='12')

    plt.subplot(1,2,2)
    ax=sb.scatterplot(x=df['l2fc'],y=df['alignment_similarity_percentage']*100,s=120, legend=False, hue =df['alignment_similarity_percentage'] ,
                   palette=sb.diverging_palette(0, 255, s=150, as_cmap=True),edgecolor='k',linewidth=0.3)

    # add gene labels
    print(df)
    df = df.reset_index(drop=True)
    print(df)
    for line in range(0,df.shape[0]):
        ax.text(df.l2fc[line], df.alignment_similarity_percentage[line]*100, 
        df.Gene[line], horizontalalignment='left', 
        size='medium', color='black', weight='semibold')
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.ylabel('Alignment Similarity %', fontsize=15, fontweight='bold')
    plt.xlabel('Log2 fold change of mean expression', fontsize = 15, fontweight='bold')
    plt.grid(False)
    # plt.axhline(50, color='black')
    plt.axvline(0, color='black', linestyle='dashed')
    plt.ylim([0,100])
    plt.savefig(path_fig + 'alignment_sim_l2fc.png')


def plot_alignmentSim_vs_optCost(x, path_fig, opt_cost_cut=0):
    sb.scatterplot(x=x['opt_alignment_cost'],y=x['alignment_similarity_percentage']*100,s=120, legend=False, hue =x['alignment_similarity_percentage'] ,
                   palette=sb.diverging_palette(0, 255, s=150, as_cmap=True),edgecolor='k',linewidth=0.3)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.ylabel('Alignment Similarity %', fontsize=15, fontweight='bold')
    plt.xlabel('Optimal alignment cost (nits)', fontsize = 15, fontweight='bold')
    plt.grid(False)
    plt.axvline(opt_cost_cut, color='black', linestyle='dashed')
    plt.tight_layout()
    plt.savefig(path_fig + 'alignment_sim.png')


def plotTimeSeriesCelltypes(al_obj, refQueryAlignerObj, cell_hue, path_fig, query, ref, gene, hue_query, hue_ref):
    plt.subplots(1,3,figsize=(15,3))
    plt.subplot(1,3,1)
    VisualUtils.plotTimeSeriesAlignment(al_obj) 
    plt.subplot(1,3,2)
    max_val = np.max([np.max(np.asarray(refQueryAlignerObj.ref_mat[al_obj.gene])), np.max(np.asarray(refQueryAlignerObj.query_mat[al_obj.gene]))])
    min_val = np.min([np.min(np.asarray(refQueryAlignerObj.ref_mat[al_obj.gene])), np.min(np.asarray(refQueryAlignerObj.query_mat[al_obj.gene]))]) 
    g = sb.scatterplot(x=refQueryAlignerObj.query_time, y=np.asarray(refQueryAlignerObj.query_mat[al_obj.gene]), alpha=0.7, hue=hue_query, legend=False,linewidth=0.3, s=20, palette=JOINT_CMAP)  
    plt.title(f'Query: {query}')
    plt.ylim([min_val-0.5,max_val+0.5])
    plt.subplot(1,3,3)
    g = sb.scatterplot(x=refQueryAlignerObj.ref_time, y=np.asarray(refQueryAlignerObj.ref_mat[al_obj.gene]), hue=hue_ref, alpha=0.7, legend=False,linewidth=0.3,s=20, palette=JOINT_CMAP)
    plt.title(f'Reference: {ref}')
    plt.ylim([min_val-0.5,max_val+0.5])
    plt.savefig(path_fig + gene +'pseudotime.png')


def get_stat_df(aligner, path_fig, FILTERED=True):

        opt_alignment_costs = [] 
        match_percentages =[] 
        l2fc = []
        for g in aligner.gene_list:
            opt_alignment_costs.append(aligner.results_map[g].fwd_DP.opt_cost)
            match_percentages.append(aligner.results_map[g].match_percentage/100)
            rgex = np.asarray(aligner.ref_mat.loc[:,g])
            qgex = np.asarray(aligner.query_mat.loc[:,g])        
            l2fc.append(np.log2(np.mean(rgex)/np.mean(qgex))) 

        df = pd.DataFrame([aligner.gene_list, match_percentages, opt_alignment_costs, l2fc]).transpose()
        df.columns = ['Gene','alignment_similarity_percentage', 'opt_alignment_cost','l2fc']
        df.set_index('Gene')
        df['color'] = np.repeat('green',df.shape[0])
        df.loc[df['alignment_similarity_percentage']<=0.5,'color'] = 'red'
        df['abs_l2fc'] = np.abs(df['l2fc']) 
        print(df)
        if FILTERED:
            df = df.sort_values(['alignment_similarity_percentage','abs_l2fc'],ascending=[True, False])

        if FILTERED:
            # remove l2fc with inf values
            df = df[~df.isin([np.inf, -np.inf]).any(axis=1)]
            df = df[df['alignment_similarity_percentage']!=0.0]
        
        plt.subplots(1,2,figsize=(16,6))
        plt.subplot(1,2,1)
        print('mean matched percentage: ')
        print(round(np.mean(df['alignment_similarity_percentage']),4)*100,'%' )
      
        temp = np.asarray([x*100 for x in df['alignment_similarity_percentage']]) 
        p = sb.kdeplot(temp, fill=True,label='Alignment Similarity %')
        plt.xlabel('Alignment Similarity Percentage')
        plt.xlim([0,100])
        p.set_yticklabels(p.get_yticks(), size = 12)
        p.set_xticklabels(p.get_xticks(), size = 12)
        plt.xlabel('Alignment similarity percentage', fontsize='12')
        plt.ylabel('Density', fontsize='12')

        plt.subplot(1,2,2)
        plot_alignmentSim_vs_optCost(df, path_fig)
        
        return df


def plot_pt(adata, root_cells, path_fig, s=10, cmap="rainbow"):
    embedding = adata.obsm['X_umap']

    pseudotime_list = [f"Pseudotime_{i}" for i in root_cells.keys()] + ["Pseudotime"]
    print("Plot pseudotime")
    print(pseudotime_list)
    for pseudotime in pseudotime_list:

        color = adata.obs[pseudotime]
        idx = np.where(~color.isna())[0]

        fig = plt.figure()
        plt.scatter(embedding[idx, 0], embedding[idx, 1], c=color[idx], s=s, cmap=cmap)
        plt.title(pseudotime)
        plt.axis("off")
        plt.savefig(path_fig + pseudotime + 'pseudotime.png', dpi=300)


def plot_root(adata, root_cells, path_fig, s=10):
    embedding = adata.obsm['X_umap']

    print("Plot root cells")
    for lineage, root_cell in root_cells.items():

        color = ["#EC7063"]
        root_cell_idx = np.where(adata.obs.index == root_cell)[0]

        fig = plt.figure()
        plt.scatter(embedding[:, 0], embedding[:, 1], color=color, s=s)
        plt.scatter(embedding[root_cell_idx, 0], embedding[root_cell_idx, 1], color="black", s=s*10, label="root_cell")
        plt.title(lineage)
        plt.axis("off")
        plt.savefig(path_fig + lineage + '_root.png', dpi=300)


# diffusion pseudotime
def pseudotime(adata, path_fig, root_cells = None, ALL=False):
    """
    Return diffusion pseudotime
    
    :param adata: Annotated scrna-seq data matrix
    :param path_fig: Save figures to this path.
    :param root_cells: optional specify the root cell for pseudotime calculation
    :param ALL: If True all celltypes will be used else only erythroid celltypes.
    """
    # Select root cells through scoring HSC marker genes
    marker_genes = ["CD34", "PROM1", "CRHBP", "NPR3", "MEIS1", "THY1", "SPINK2"]
    sc.tl.score_genes(adata, marker_genes)

    print("Pseudotime calculation\n")
    # Instantiate pseudotime object using anndata object.
    pt = Pseudotime_calculator(adata=adata,
                            obsm_key="X_umap", # Dimensional reduction data name
                            cluster_column_name=CELL_CLUSTER # Clustering data name
                            )

    # change lineage dict in case we calculate PT across entire dataset
    if ALL:
        lineage_dictionary = LINEAGE_DICT
    else:
        lineage_dictionary = {"Lineage_Eryth": CT_ERY}

    # Input lineage information into pseudotime object
    pt.set_lineage(lineage_dictionary=lineage_dictionary)

    # find and set root cells
    if ALL:
        adata_tmp = adata[adata.obs[CELL_CLUSTER].isin(['HSCs & MPPs'])]
        if adata_tmp.n_obs > 0:
            print(f"Index of cells with the highest scores {np.argmax(adata_tmp.obs['score'])}")
            hsc_root = adata_tmp.obs.index[np.argmax(adata_tmp.obs['score'])]
            b_t_nk_root = hsc_root
            print(f"The hsc root is {hsc_root}.")
        else: 
            # if there are no stem cells in the set use EMP and LMP as roots for the lineage
            adata_tmp = adata[adata.obs[CELL_CLUSTER].isin(['EMP'])]
            print(f"Index of EMP cells with the highest scores {np.argmax(adata_tmp.obs['score'])}")
            hsc_root = adata_tmp.obs.index[np.argmax(adata_tmp.obs['score'])]
            print(f"The EMP root is {hsc_root}.")
            adata_tmp = adata[adata.obs[CELL_CLUSTER].isin(['LMP'])]
            print(f"Index of cells with the highest scores {np.argmax(adata_tmp.obs['score'])}")
            b_t_nk_root = adata_tmp.obs.index[np.argmax(adata_tmp.obs['score'])]
            print(f"The B/ T Nk root is {b_t_nk_root}.")
        
        root_cells = {"Lineage_Eryth": hsc_root,
                            "Lineage_Thrombo": hsc_root,
                            "Lineage_B": b_t_nk_root,
                            "Lineage_T_Nk": b_t_nk_root,
                            "Lineage_Dendritic": hsc_root,
                            "Lineage_Myelo": hsc_root,
                            }
    else:
        root_cells = {"Lineage_Eryth": hsc_root}

    print("Set root cells")
    pt.set_root_cells(root_cells=root_cells)

    # Check diffusion map data.
    if "X_diffmap" not in pt.adata.obsm:
        print("Calculate diffusion map")
        # Calculate diffusion map if anndata object does not have diffusion map data.
        sc.pp.neighbors(pt.adata, n_neighbors=30)
        sc.tl.diffmap(pt.adata)

    # Calculate pseudotime
    try:
        pt.get_pseudotime_per_each_lineage()
    except ValueError as ve:
        print("Error message:")
        print(ve)

        # plot igraph
        nbrs =sc.Neighbors(pt.adata)
        g = nbrs.to_igraph()
        fig, ax = plt.subplots()
        ig.plot(g, target=ax)
        plt.savefig(path_fig +'neighbors_igraph.png')
        return -1

    plot_pt(pt.adata, root_cells, path_fig)
    plot_root(pt.adata, root_cells, path_fig + f"{ROOT_ID}_")
    adata.obs['time'] = pt.adata.obs[["Pseudotime"]]
    return adata


# Plot pseudotime trajectory alignment of three groups at once
def plot_three_alignments(aligner_1, aligner_2, gene, path_fig, ref, query_1, query_2):
    """
    Returns plot showing alignment for three groups. 
    Two queries are aligned against the same reference.

    :param aligner_1: Alignement object of ref and query_1
    :param aligner_2: Alignement object of ref and query_2
    :param gene: String label for gene of the alignment
    :param path_fig: Save figures to path_fig.
    :param ref: String label for ref
    :param query_1: String label for query_1
    :param query_2: String label for query_2
    """
    legend_elements = [Line2D([0], [0], marker='o', color='forestgreen', label=ref,
                            markerfacecolor='forestgreen', markersize=5),
                    Line2D([0], [0], marker='o', color='midnightblue', label=query_1,
                            markerfacecolor='midnightblue', markersize=5),
                    Line2D([0], [0], marker='o', color='lightblue', label=query_2,
                            markerfacecolor='lightblue', markersize=5),
                    ]
    
    fig, ax = plt.subplots()
    al_obj_1 = aligner_1.results_map[gene]
    al_obj_2 = aligner_2.results_map[gene]
    sb.scatterplot(x=al_obj_1.S.X, y=al_obj_1.S.Y, color = 'forestgreen' ,alpha=0.05, legend=False) #Ref
    sb.scatterplot(x=al_obj_1.T.X, y=al_obj_1.T.Y, color = 'midnightblue' ,alpha=0.05, legend=False) #Query 1
    sb.scatterplot(x=al_obj_2.T.X, y=al_obj_2.T.Y, color = 'lightblue' ,alpha=0.05, legend=False) #Query 2
     
    al_obj_1.S.plot_mean_trend(color='forestgreen') #Ref
    al_obj_1.T.plot_mean_trend(color='midnightblue') #Query 1
    al_obj_2.T.plot_mean_trend(color='lightblue') #Query 2
    plt.title(gene, fontsize = 20)
    plt.xlabel('Pseudotime', fontsize = 16)
    plt.ylabel('Gene expression (log-normalized)', fontsize = 16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)

    ax.legend(handles=legend_elements, fontsize = 16,loc='center left', bbox_to_anchor=(1, 0.5))

    sub_path_fig = f"{path_fig}/_{ref}_{query_1}_{query_2}/"  
    Path(sub_path_fig).mkdir(parents=True, exist_ok=True)

    plt.savefig(sub_path_fig + gene + '.png', dpi=300, bbox_inches='tight')
    return fig, ax


def plot_distances(df, path_fig, ref, query):
    plt.rcParams.update({'font.size': 18})
    fig = plt.figure(figsize =(10, 7))
    for f_dist,y in zip(df["frechet_distance"], range(df.shape[0])):
        plt.plot(f_dist, y, 'g*')
    plt.yticks(range(df.shape[0]),list(df.index))
    plt.xlim([0, np.max(df['frechet_distance'])+ 0.2] )
    plt.title(f"Frechet distance {ref} vs {query}")
    plt.savefig(path_fig +'frechet_dist_barplots.png')

    plt.rcParams.update({'font.size': 18})
    fig = plt.figure(figsize =(10, 7))
    for dist,y in zip(df["auc_dist"], range(df.shape[0])):
        plt.plot(dist, y, 'g*')
    plt.yticks(range(df.shape[0]),list(df.index))
    plt.xlim([0, np.max(df['auc_dist'])+ 0.2] )
    plt.title(f"AUC distance {ref} vs {query}")
    plt.savefig(path_fig +'auc_dist_barplots.png')

    plt.rcParams.update({'font.size': 18})
    fig = plt.figure(figsize =(10, 7))
    for dist,y in zip(df["mean_abs_dist"], range(df.shape[0])):
        # plt.plot((lower,mean, upper),(y,y,y), 'b.-')
        plt.plot(dist, y, 'g*')
    plt.yticks(range(df.shape[0]),list(df.index))
    plt.xlim([0, np.max(df['mean_abs_dist'])+ 0.2] )
    plt.title(f"Mean absolute distance {ref} vs {query}")
    plt.savefig(path_fig +'mean_abs_dist_barplots.png')

    return


def visualize_CI(path_fig, file, metric, genes, title_name =""):
    df = pd.read_csv(file, index_col=0) 
    # sort genes
    df = df.loc[genes]
    df = df.reindex(genes[::-1])
    plt.rcParams.update({'font.size': 18})
    fig = plt.figure(figsize =(10, 7))
    for lower, upper, mean, orig, y in zip(df['CI_lower'], df['CI_upper'], df['mean'], df[f"original_{metric}"], range(df.shape[0])):
        plt.plot((lower,mean, upper),(y,y,y), 'b.-')
        plt.plot(orig, y, 'g*')
    plt.yticks(range(df.shape[0]),list(df.index))
    plt.xlim([0, max(df['CI_upper'])+0.2])
    plt.xlabel("Expression dissimilarity")
    if title_name == "":
        plt.title(metric)
    else:
        plt.title(title_name)
    plt.savefig(path_fig +f'{metric}_CI_barplots.png')
    

def ensure_unit_interval(time_array):
    """
    Return min-max normalized pseudotime if not normalized yet
    
    :param time_array: Pseudotime
    """
    arr = np.asarray(time_array)
    if not (np.isclose(arr.min(), 0.0) and np.isclose(arr.max(), 1.0)):
        arr = TimeSeriesPreprocessor.Utils.minmax_normalise(arr)
    return arr


def threewise_alignment(adata_ref, adata_query_1, adata_query_2, genes, ref, query_1, query_2, path_fig):
    """
    Perform pseudotime gene expression alignment of three groups.
    
    :param adata_ref: Annotated data matrix for reference
    :param adata_query_1: Annotated data matrix for query_1
    :param adata_query_2: Annotated data matrix for query_2
    :param genes: Alignment is performed for these genes
    :param ref: String label for ref
    :param query_1: String label for query_1
    :param query_2: String label for query_2
    :param path_fig: Save figures to path_fig
    """
    print("G2G alignment with threewise plot\n")

    # if not yet normalized min-max normalize pseudotime
    adata_ref.obs['time'] = ensure_unit_interval(adata_ref.obs['time'])
    adata_query_1.obs['time'] = ensure_unit_interval(adata_query_1.obs['time'])
    adata_query_2.obs['time'] = ensure_unit_interval(adata_query_2.obs['time'])

    # determine number of bins / interpolation points for alignment
    x = np.asarray(adata_ref.obs.time)
    optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
    optb.fit(x, x)
    opt_splits_ref = len(optb.splits)
    print(f"the optimal number of splits for {ref} is {len(optb.splits)}")

    x = np.asarray(adata_query_1.obs.time)
    optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
    optb.fit(x, x)
    opt_splits_query_1 = len(optb.splits)
    print(f"the optimal number of splits for {query_1} is {len(optb.splits)}")

    x = np.asarray(adata_query_2.obs.time)
    optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
    optb.fit(x, x)
    opt_splits_query_2 = len(optb.splits)
    print(f"the optimal number of splits for {query_2} is {len(optb.splits)}")

    n_bins = max(opt_splits_ref, opt_splits_query_1, opt_splits_query_2)
    print(f"The optimal number of splits selected for the alignment is {n_bins}")

    print(f'Aligning {ref} against {query_1}')
    aligner_1 = Main.RefQueryAligner(adata_ref, adata_query_1, genes, n_bins)
    aligner_1.WEIGHT_BY_CELL_DENSITY = True
    aligner_1.WINDOW_SIZE=0.1
    aligner_1.align_all_pairs()

    print(f'Aligning {ref} against {query_1}')
    aligner_2 = Main.RefQueryAligner(adata_ref, adata_query_2, genes, n_bins)
    aligner_2.WEIGHT_BY_CELL_DENSITY = True
    aligner_2.WINDOW_SIZE=0.1
    aligner_2.align_all_pairs()

    for g in genes:
        plot_three_alignments(aligner_1, aligner_2, g, path_fig, ref, query_1, query_2)

    return 


# use hd base to normalize gene expression according to healthy donors
def g2g_alignment(adata_ref, adata_query, hd_base, path_fig, genes, query, ref, plot_genes = False, FILTERED=True, DIST = True, EQ =True):
    """
    Pseudotime gene expression alignment based on genes2genes.
    
    :param adata_ref: Annotated matrix for reference
    :param adata_query: Annotated matrix for query
    :param hd_base: Normalize gene expression against healthy donors
    :param path_fig: Figures are saved at path_fig
    :param genes: Alignment is performed for genes
    :param query: String label for query
    :param ref: String label for ref
    :param plot_genes: If True plot genewise alignment
    :param FILTERED: If True filter out infiniy values
    :param F_DIST: If True calculate distance in alignment
    :param EQ: If True set pseudotime to same length
    """
    print("G2G alignment\n")

    # if not yet normalized min-max normalize pseudotime
    adata_ref.obs['time'] = ensure_unit_interval(adata_ref.obs['time'])
    adata_query.obs['time'] = ensure_unit_interval(adata_query.obs['time'])

    # reduce data size by only keeping info on genes that we need
    hd_base = hd_base[:, genes]

    # prep average gene expression needed for normalizing
    hd_avg_gene_exp =  pd.DataFrame(hd_base.X.todense()) 
    hd_avg_gene_exp.columns = hd_base.var_names
    hd_avg_gene_exp = hd_avg_gene_exp.set_index(hd_base.obs_names)
    hd_avg_gene_exp =hd_avg_gene_exp.mean(axis=0)
    print(f"Mean genewise expression of HD: {hd_avg_gene_exp}")

    if EQ:
        print("Equalize pseudotime")
        time_min = min(max(adata_ref.obs['time']), max(adata_query.obs['time']))
        print(f"Minimum maximum time {time_min}")

        adata_ref = adata_ref[adata_ref.obs['time'] <= time_min]
        adata_query = adata_query[adata_query.obs['time'] <= time_min]

    # determine number of bins / interpolation points for alignment
    x = np.asarray(adata_ref.obs.time)
    optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
    optb.fit(x, x)
    opt_splits_ref = len(optb.splits)
    print(f"the optimal number of splits for {ref} is {len(optb.splits)}")

    x = np.asarray(adata_query.obs.time)
    optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
    optb.fit(x, x)
    opt_splits_query = len(optb.splits)
    print(f"the optimal number of splits for {query} is {len(optb.splits)}")

    n_bins = max(opt_splits_ref, opt_splits_query)
    print(f"The optimal number of splits selected for the alignment is {n_bins}")

    # Aligning all genes in the given list at once
    aligner = Main.RefQueryAligner(adata_ref, adata_query, genes, n_bins)
    aligner.WEIGHT_BY_CELL_DENSITY = True
    aligner.WINDOW_SIZE=0.1
    aligner.align_all_pairs()

    if plot_genes:
        # create subfolder for the timeseries
        Path(path_fig).mkdir(parents=True, exist_ok=True)
        sub_path_fig = path_fig + "/" #requires at the moment that the new directory name is already in the pathfig
        for g in genes:
            plotTimeSeriesCelltypes(aligner.results_map[g], aligner, 
                                    cell_hue= CELL_CLUSTER, 
                                    gene = g, 
                                    query = query, 
                                    ref = ref, 
                                    path_fig = sub_path_fig,
                                    hue_query = adata_query.obs[CELL_CLUSTER],
                                    hue_ref = adata_ref.obs[CELL_CLUSTER])
            
    
    frechet_dist, auc_dist, mean_abs_dist = [], [], []
    df = None
    if DIST:
        print("Calculate distances")
        for g in genes:
            min_len = min(len(aligner.results_map[g].S.mean_trend), len(aligner.results_map[g].T.mean_trend))
            ref_mean_trend = aligner.results_map[g].S.mean_trend[:min_len]
            query_mean_trend = aligner.results_map[g].T.mean_trend[:min_len]
            ref_time = aligner.results_map[g].S.time_points[:min_len]
            query_time = aligner.results_map[g].T.time_points[:min_len]

            # normalize gene expression by dividing it by the healthy baseline
            avg_hd = hd_avg_gene_exp[g]
            ref_mean_trend = ref_mean_trend /avg_hd
            query_mean_trend = query_mean_trend /avg_hd
            

            ref_line = LineString(list(zip(ref_time, ref_mean_trend)))
            query_line = LineString(list(zip(query_time, query_mean_trend)))
            frechet_dist.append(shapely.frechet_distance(ref_line, query_line))

            #auc distance = abs(trapz_auc_ref -trapz_auc_query)
            ref_auc = np.trapz(ref_mean_trend)
            query_auc = np.trapz(query_mean_trend)
            auc_dist.append(abs(ref_auc - query_auc))

            # mean abs dist = abs(y_1-y_2)
            mean_abs_dist.append(np.mean(abs(ref_mean_trend - query_mean_trend)))


        df = pd.DataFrame([genes, frechet_dist, auc_dist, mean_abs_dist]).transpose()
        df.columns = ['Gene','frechet_distance', 'auc_dist', 'mean_abs_dist']
        df.set_index('Gene', inplace=True) 
        print(df)
        df.to_csv(path_fig + "distances_df.csv")

        plot_distances(df, path_fig, ref, query)
        
    else: # alignment similarity case
        df = get_stat_df(aligner, path_fig, FILTERED) 
        
        plot_alignmentSim_vs_l2fc(df, path_fig)
        df.to_csv(path_fig + "stat_df.csv")


    return df


def bootstrapping(adata_ref, adata_query, hd_base, n_iter, path_fig, genes, query, ref):
    """
    Perform bootstrapping of pseudotime alignment
    
    :param adata_ref: Annotated matrix for reference
    :param adata_query: Annotated matrix for query
    :param hd_base: Normalize gene expression against healthy donors
    :param n_iter: Resample n_iter times
    :param path_fig: Save figures at path_fig
    :param genes: Alignment is performed for genes
    :param query: String label for query
    :param ref: String label for ref
    """
    fre_dist = pd.DataFrame(np.zeros((len(genes), n_iter)))
    auc_dist = pd.DataFrame(np.zeros((len(genes), n_iter)))
    am_dist = pd.DataFrame(np.zeros((len(genes), n_iter)))
    fre_dist.index, auc_dist.index, am_dist.index = genes, genes, genes

    # Print original cell type proportion 
    print(f"Cell types ref: {print(adata_ref.obs[CELL_CLUSTER].value_counts())}")
    print(f"Cell types query: {print(adata_query.obs[CELL_CLUSTER].value_counts())}")

    i = 0
    while (i<n_iter):
        print(f"Resampling round: {i}")
        np.random.seed()

        # resample with replacement for each group, resample from counts
        # sampling mask, fixing that the root cell is in the sample
        ref_mask = np.random.randint(adata_ref.n_obs, size=adata_ref.n_obs)
        query_mask = np.random.randint(adata_query.n_obs, size=adata_query.n_obs)

        adata_ref_rwr = adata_ref[ref_mask, :]
        adata_query_rwr = adata_query[query_mask, :]

        # reset original names for pseudotime calculations
        adata_ref_rwr.obs_names = adata_ref.obs_names
        adata_query_rwr.obs_names = adata_query.obs_names

        if adata_ref_rwr is None or adata_query_rwr is None:
            print("Failed pseudotime => Redo calculation")
            continue
        else: 
            # calculate alignment
            df = g2g_alignment(adata_ref = adata_ref_rwr, 
                               adata_query = adata_query_rwr, 
                               hd_base= hd_base,
                               path_fig = f"{path_fig}bst_", 
                               genes = genes, 
                               query = query, 
                               ref = ref, 
                               FILTERED=False, 
                               plot_genes = False)
            fre_dist.iloc[:, i] = df['frechet_distance'].values
            auc_dist.iloc[:, i] = df['auc_dist'].values
            am_dist.iloc[:, i] = df['mean_abs_dist'].values
            i += 1

    # original result
    df = g2g_alignment(adata_ref = adata_ref, adata_query = adata_query, hd_base= hd_base, 
                       path_fig = f"{path_fig}", genes = genes, query = query, ref = ref, FILTERED=False)

    for met,dist_df in zip(['frechet_distance', 'auc_dist', 'mean_abs_dist'], [fre_dist, auc_dist, am_dist]):    
        # calculate mean / median of alignment similarity 
        means = np.mean(dist_df, axis = 1)
        confidence_interval = np.percentile(dist_df, [2.5, 97.5], axis = 1)
        dist_df[f"original_{met}"] = df[met].values
        dist_df["mean"] = means
        dist_df["CI_lower"] = confidence_interval[0]
        dist_df["CI_upper"] = confidence_interval[1]
        dist_df.to_csv(path_fig + f"{met}_bootstrapped_samples.csv")
        visualize_CI(path_fig, path_fig + f"{met}_bootstrapped_samples.csv", met, genes)


def store_pseudotime(path_fig, file_name, PT_cases, selection_col, results_path, SAVE = False, ALL_CELLS = True):
    """
    Save pseudotime in adata object
    
    :param path_fig: Save figures at path_fig
    :param file_name: Read in adata from file_name
    :param PT_cases: Store pseudotime for each case in PT_cases
    :param selection_col: Observation colum for PT case
    :param results_path: If SAVE =True save adata at results_path
    :param SAVE: If True save adata object
    :param ALL_CELLS: If True use all celltypes else only consider erythroid cells
    """
    # set directory and read in adata object, which has been normalized and log-transformed
    sc.settings.figdir = path_fig
    Path(path_fig).mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(file_name)

    # save pseudotime => PT stored in adata.obs['time']
    # iterate all comparisons, calculcate pseudotime and after save adata object
    for case, col in zip(PT_cases, selection_col):
        print(f"Pseudotime case: {case} in column {col}\n")
        print(adata.obs[col])

        if case == "ED MDS/AML-clus1":
            adata_case = adata[adata.obs[col].isin(["ED MDS/AML"])]
            adata_case = adata_case[adata_case.obs['leiden_EMP_0.06']!= '1']
        else:
            adata_case = adata[adata.obs[col].isin([case])] 
        print(f"There are {adata_case.n_obs} cells that are considered to be {case}.")
        adata_case = pseudotime(adata_case, path_fig = f"{path_fig}case_", ALL = ALL_CELLS)
        print(adata_case)

        # store PT for each condition 
        if case == "ED MDS/AML-clus1" or case == "ED MDS/AML":
            case = case.replace("/", " wo ") # adata does not like backslash in obs names
        adata.obs[f'PT_{case}'] = np.nan
        adata.obs.loc[adata_case.obs_names, f'PT_{case}'] = adata_case.obs['time']
    print(adata)    
    
    if SAVE:
        adata.write(results_path + 'adata_pt.h5ad')



def pt_comparisons(path_base, file_name, refs, queries, selection_col, celltypes_to_keep, marker_genes, CT):
    """
    Perform pseutotime gene expression alignment for given comparisons
    
    :param path_base: Save figures at path_base
    :param file_name: Read in adata from file_name
    :param refs: References for comparisons
    :param queries: Queries for comparisons
    :param selection_col: adata.obs column to find refs and queries
    :param celltypes_to_keep: Use celltypes for comparisons
    :param marker_genes: Alignment is performed for marker_genes
    :param CT: Label for celltypes that we are investigating
    """
    assert len(queries)==len(refs), "There are different number of cases in refs and queries"

    adata = sc.read_h5ad(file_name)
    adata = adata[adata.obs[CELL_CLUSTER].isin(celltypes_to_keep)]
    
    # if one combined pseudotime for all these, uncomment the following lines
    # case = ["HD", "ED BMF", "ED MDS", "ED BM normal","SDS BMF", "SDS MDS", "ED MDS/AML"]
    # col = 'condition_diagnosis'
    # adata = adata[adata.obs[col].isin(case)]
    # print(f"Normalize PT_avg_combined")
    # adata.obs['PT_avg_combined'] = ensure_unit_interval( adata.obs['PT_avg_combined'])

    # iterate through reference and query sets and calculate alignments 
    for ref, query, col in zip(refs, queries, selection_col):
        print(f"Alignment for {ref} versus {query} in Column {col}")

        # prepare name of EML/MDS-clus-1
        query_dir = query.replace(" ", "_")
        query_dir = query_dir.replace("/", " wo ")
        # set figure directories
        path_fig = path_base + "average_PT/" + f'ref_{ref.replace(" ", "_")}_query_{query_dir}/{CT}/' 
        print(path_fig)
        sc.settings.figdir = path_fig
        Path(path_fig).mkdir(parents=True, exist_ok=True)

        #special case HD HCA, stored differently in the data, should only exist in the references
        if ref == "HD_HCA":
            adata_ref = adata[adata.obs['individual'].str.startswith('HCA_')]
        else:
            adata_ref = adata[adata.obs[col].isin([ref])] 

        adata_ref.obs['time'] = adata_ref.obs[f'PT_avg_{ref.replace("/", "wo") }'] # / were removed in the obs names
        # adata_ref.obs['time'] = adata_ref.obs[f'PT_avg_combined'] # combined PT for all

        if query == "ED MDS/AML-clus1":
            adata_query = adata[adata.obs[col].isin(["ED MDS/AML"])]
            adata_query = adata_query[adata_query.obs['leiden_EMP_0.06']!= '1']
        else:
            adata_query = adata[adata.obs[col].isin([query])] 
        adata_query.obs['time'] = adata_query.obs[f'PT_avg_{query.replace("/", "wo") }'] # / were removed in the obs names
        # adata_query.obs['time'] = adata_query.obs[f'PT_avg_combined'] # combined PT for all

        g2g_alignment(adata_ref = adata_ref,
            adata_query = adata_query,
            hd_base= adata[adata.obs[col].isin(["HD"])],
            path_fig = f"{path_fig}marker_", 
            genes = marker_genes,
            query = query,
            ref = ref, 
            plot_genes = False)  
        
    return


def bootstrapping_pt_comp(path_base, file_name, refs, queries, selection_col, celltypes_to_keep, marker_genes):
    """
    Setup bootstrapping of pseudotime alignemnts for given comparisons

    :param path_base: Save figures at path_base
    :param file_name: Read in adata from file_name
    :param refs: References for comparisons
    :param queries: Queries for comparisons
    :param selection_col: adata.obs column to find refs and queries
    :param celltypes_to_keep: Use celltypes for comparisons
    :param marker_genes: Alignment is performed for marker_genes
    """
    assert len(queries)==len(refs), "There are different number of cases in refs and queries"
    adata = sc.read_h5ad(file_name)
    
    print(f"Marker genes: {marker_genes}")

    adata = adata[adata.obs[CELL_CLUSTER].isin(celltypes_to_keep)]
    case = ["HD", "ED BMF", "ED MDS", "ED BM normal","SDS BMF", "SDS MDS", "ED MDS/AML"]
    col = 'condition_diagnosis'
    adata = adata[adata.obs[col].isin(case)]

    print(f"Normalize PT_avg_combined")
    adata.obs['PT_avg_combined'] = ensure_unit_interval( adata.obs['PT_avg_combined'])

    # iterate through reference and query sets and calculate alignments 
    for ref, query, col in zip(refs, queries, selection_col):
        print(f"Confidence interval for alignment for {ref} versus {query}")
        # prepare name of EML/MDS-clus-1
        query_dir = query.replace(" ", "_")
        query_dir = query_dir.replace("/", " wo ")
        ref_dir = ref.replace(" ", "_")
        ref_dir = ref_dir.replace("/", " wo ")
        # set figure directories
        path_fig = path_base + "average_PT/" + f'ref_{ref_dir}_query_{query_dir}/' 
        print(path_fig)
        sc.settings.figdir = path_fig
        Path(path_fig).mkdir(parents=True, exist_ok=True)

        #special case HD HCA, stored differently in the data, should only exist in the references
        if ref == "HD_HCA":
            adata_ref = adata[adata.obs['individual'].str.startswith('HCA_')]
        else:
            adata_ref = adata[adata.obs[col].isin([ref])] 

        adata_ref.obs['time'] = adata_ref.obs[f'PT_avg_combined'] # combined PT for all
        print(adata_ref.obs['time'])

        if query == "ED MDS/AML-clus1":
            adata_query = adata[adata.obs[col].isin(["ED MDS/AML"])]
            adata_query = adata_query[adata_query.obs['leiden_EMP_0.06']!= '1']
        else:
            adata_query = adata[adata.obs[col].isin([query])] 
        adata_query.obs['time'] = adata_query.obs[f'PT_avg_combined'] # combined PT for all

        bootstrapping(adata_ref, 
                  adata_query, 
                  hd_base= adata[adata.obs['condition_diagnosis'].isin(["HD"])],
                  n_iter = 100, 
                  path_fig = f"{path_fig}", 
                  genes = marker_genes,
                  ref = ref,
                  query = query)    
        for metric in ['auc_dist', 'mean_abs_dist']:
            if ref=="HD":
                ref = "HC"
            title_name = f"{query} vs {ref}"
            visualize_CI(f"{path_fig}marker", f"{path_fig}marker_" + f"{metric}_bootstrapped_samplesmerged.csv", metric, marker_genes, title_name)

    return


def pseudotime_bootstrapping(path_fig, case, selection_col, file_name, pt_file_name):
    """
    Perform bootstrapping of pseudotime
    
    :param path_fig: Save figures at path_fig
    :param case: Case/s for which pseudotime is calculated
    :param col: adata.obs column to find cases
    :param file_name: Read in adata from file_name
    :pt_file_name: Save bootstrapped pseudotimes under this pt_file_name
    """
    adata_case = sc.read_h5ad(file_name)
    
    if case == "ED MDS/AML-clus1":
            adata_case = adata_case[adata_case.obs[selection_col].isin(["ED MDS/AML"])]
            adata_case = adata_case[adata_case.obs['leiden_EMP_0.06']!= '1']
    elif len(case)>1: # PT for combined cases
        adata_case = adata_case[adata_case.obs[selection_col].isin(case)] 
    else:
        adata_case = adata_case[adata_case.obs[selection_col].isin([case])] 
    print(f"There are {adata_case.n_obs} cells that are considered to be {case}.")

    # set root cells for pseudotime
    num_hsc = len(adata_case[adata_case.obs[CELL_CLUSTER].isin(['HSCs & MPPs'])]) # number of stem cells
    EMP_FLAG = False

    if num_hsc == 0: # if there are no stem cells iterate through emps 
        num_hsc = len(adata_case[adata_case.obs[CELL_CLUSTER].isin(['EMP'])])
        EMP_FLAG = True

    # separate case for HD and HD TP53 WT or if we look at combined cases
    # randomly sample from stem cells as many more stem cells exist for this case
    if case in ["HD", "HD TP53 WT"] or len(case)>1:
        num_hsc_sample = 5 # highest number of stem cells observed in other cases

    np.random.seed()
    hsc_mask = np.random.randint(num_hsc, size=num_hsc_sample) # resample with replacement 
    print(f"Subsampling for {case} - the resampled root cell indices are {hsc_mask}")
    num_hsc = num_hsc_sample

    print(f"There exist {num_hsc} stem cells in the {case} data.")
    sys.stdout.flush()
    pt_df = pd.DataFrame(np.zeros((int(adata_case.n_obs), int(num_hsc))))
    pt_df.index = adata_case.obs_names # set cell ids as index
   
    # Instantiate pseudotime object using anndata object.
    pt = Pseudotime_calculator(adata=adata_case,
                            obsm_key="X_umap", # Dimensional reduction data name
                            cluster_column_name=CELL_CLUSTER # Clustering data name
                            )
    # Input lineage information into pseudotime object
    pt.set_lineage(lineage_dictionary=LINEAGE_DICT)

    # iterate through HSCs as root cells and calculate pseudotime
    for i in range(num_hsc):
        if not EMP_FLAG:
            hsc_root = adata_case.obs.loc[lambda x: x[CELL_CLUSTER] == "HSCs & MPPs"].index[i]
            b_t_nk_root = hsc_root
            print(f"The hsc root is {hsc_root}.")
        elif case in ["HD", "HD TP53 WT"]:
            print("HD/ HD TP53 WT case")
            sampled_i = hsc_mask[i]
            hsc_root = adata_case.obs.loc[lambda x: x[CELL_CLUSTER] == "HSCs & MPPs"].index[sampled_i]
            b_t_nk_root = hsc_root
            print(f"The hsc root is {hsc_root}.")
        elif len(case)>1:
            print("Combined cases")
            sampled_i = hsc_mask[i]
            adata_hd =  adata_case[adata_case.obs[selection_col].isin(['HD'])]
            hsc_root = adata_hd.obs.loc[lambda x: x[CELL_CLUSTER] == "HSCs & MPPs"].index[sampled_i] 
            b_t_nk_root = hsc_root
            print(f"The hsc root is {hsc_root}.")
        else: 
            # if there are no stem cells in the set use EMP and LMP as roots for the lineage
            hsc_root = adata_case.obs.loc[lambda x: x[CELL_CLUSTER] == "EMP"].index[i]
            print(f"The EMP root is {hsc_root}.")
            b_t_nk_root = adata_case.obs.loc[lambda x: x[CELL_CLUSTER] == "LMP"].index[i]
            print(f"The B/ T Nk root is {b_t_nk_root}.")


        root_cells = {"Lineage_Eryth": hsc_root,
                        "Lineage_Thrombo": hsc_root,
                        "Lineage_B": b_t_nk_root, 
                        "Lineage_T_Nk": b_t_nk_root,
                        "Lineage_Dendritic": hsc_root,
                        "Lineage_Myelo": hsc_root,
                        }

        print("Set root cells")
        pt.set_root_cells(root_cells=root_cells)

        # Check diffusion map data.
        if "X_diffmap" not in pt.adata.obsm:
            print("Calculate diffusion map")
            # Calculate diffusion map if your anndata object does not have diffusion map data.
            sc.pp.neighbors(pt.adata, n_neighbors=30)
            sc.tl.diffmap(pt.adata)

        # Calculate pseudotime
        try:
            pt.get_pseudotime_per_each_lineage()
            pt_df.iloc[:, i] = pt.adata.obs[["Pseudotime"]]
        except ValueError as ve:
            print("Error message:")
            print(ve)
            print(pt.adata.obs["Pseudotime"])
            print(pt.adata.obsm["X_diffmap"])
            print(pt.adata.uns["neighbors"])

            # plot igraph
            nbrs =sc.Neighbors(pt.adata)
            g = nbrs.to_igraph()
            fig, ax = plt.subplots()
            ig.plot(g, target=ax)
            plt.savefig(path_fig +f'neighbors_igraph_root{i}.png')
            continue
            
    pt_df.to_csv(pt_file_name)
    return 


def save_average_PT(cases, cols, file_name, pt_file_name, result_path, NORM=False):
    """
    Save average pseudotime from bootstrapped pseudotimes.
    
    :param cases: Save average pseudotime for these cases
    :param cols: adata.obs columns to find cases
    :param file_name: Read in adata from file_name
    :param pt_file_name: Read in bootstrapped pseudotime results from pt_file_name
    :param result_path: Save adata at result_path
    :param NORM: If True min-max normalize pseudotime
    """
    adata = sc.read_h5ad(file_name)
    
    for case, col in zip(cases, cols):
        if case == "ED MDS/AML-clus1":
            adata_case = adata[adata.obs[col].isin(["ED MDS/AML"])]
            adata_case = adata_case[adata_case.obs['leiden_EMP_0.06']!= '1']
        else:
            adata_case = adata[adata.obs[col].isin(case)] 
        print(f"There are {adata_case.n_obs} cells that are considered to be {case}.")

        pt_df = pd.read_csv(pt_file_name, index_col=0, header=0) 
        adata_case.obs[f"PT_combined_average"] = pt_df.mean(axis=1)

        if len(case)>1:
            case = "combined"

        if NORM:
            print(f"Normalize 'PT_{case}_average'")
            adata_case.obs[f"PT_{case}_average"] = ensure_unit_interval(adata_case.obs[f"PT_{case}_average"])
    
        ci = pt_df.apply(lambda x: stats.t.interval(0.95, len(x)-1, loc=np.mean(x), scale=stats.sem(x)), axis=1)
        ci = pd.DataFrame(ci.tolist(), index=ci.index, columns=['lb', 'ub'])
        adata_case.obs[f"PT_{case}_lb"] = ci['lb']
        print(adata_case.obs[f"PT_{case}_lb"])
        adata_case.obs[f"PT_{case}_ub"] = ci['ub']
        print(adata_case.obs[f"PT_{case}_ub"])
        sys.stdout.flush()

        adata.obs[f'PT_avg_{case}'] = np.nan
        adata.obs.loc[adata_case.obs_names, f'PT_avg_{case}'] = adata_case.obs[f'PT_{case}_average']
        print(adata.obs[f'PT_avg_{case}'])
    
    print(adata)
    adata.write(result_path)


def trio_comparison(path_base, file_name, celltypes_to_keep, marker_genes):
    """
    Perfom alignments for three groups at once across different comparisons
    
    :param path_base: Save figures at path_base
    :param file_name: Read in adata from file_name
    :param celltypes_to_keep: Use celltypes for comparisons
    :param marker_genes: Alignment is performed for marker_genes
    """
    adata = sc.read_h5ad(file_name)

    refs = ["HD", "HD TP53 WT", "HD TP53 WT"]
    queries_ed = ["ED BMF", "ED BMF TP53 WT", "ED TP53 WT"]
    queries_sds = ["SDS BMF", "SDS BMF TP53 WT", "SDS TP53 WT"]

    selection_col = ['condition_diagnosis', 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt_BMFMDS', 'tp53_genotype_wtmut_noTP53clone_wt_HD_wt', ]
    adata = adata[adata.obs[CELL_CLUSTER].isin(celltypes_to_keep)]

    for ref, query_ed, query_sds, col in zip(refs, queries_ed, queries_sds, selection_col):
        print(f"Trio plotting for {ref}, {query_ed} and {query_sds}")

        # set figure directories
        path_fig = path_base + "trio_alignments/"
        print(path_fig)
        sc.settings.figdir = path_fig
        Path(path_fig).mkdir(parents=True, exist_ok=True)

        adata_ref = adata[adata.obs[col].isin([ref])]
        adata_ref.obs['time'] = adata_ref.obs[f'PT_avg_{ref.replace("/", "wo") }'] # / were removed in the obs names

        adata_query_ed = adata[adata.obs[col].isin([query_ed])]
        adata_query_ed.obs['time'] = adata_query_ed.obs[f'PT_avg_{query_ed.replace("/", "wo") }'] # / were removed in the obs names

        adata_query_sds = adata[adata.obs[col].isin([query_sds])]
        adata_query_sds.obs['time'] = adata_query_sds.obs[f'PT_avg_{query_sds.replace("/", "wo") }'] # / were removed in the obs names

        threewise_alignment(adata_ref,
                            adata_query_1 = adata_query_ed,
                            adata_query_2 = adata_query_sds,
                            genes = marker_genes,
                            ref = ref,
                            query_1 = query_ed,
                            query_2 = query_sds,
                            path_fig = path_base + "trio_alignments/" )


    return


def sixwise_comp(path_base, file_name, EQ, marker_genes, celltypes):
    """
    Plotting pseudotime gene expression comparison across the six different conditions
    
    :param path_base: Save figures at path_base
    :param file_name: Read in adata from file_name
    :param EQ: If True set pseudotime to same length
    :param marker_genes: Alignment is performed for marker_genes
    :param celltypes: Use celltypes for comparisons
    """
    adata = sc.read_h5ad(file_name)
    selection_col = 'condition_diagnosis'
    
    # need to add GATA2/GATA1 obs to adata
    if 'GATA2/GATA1' in marker_genes:
        gata2_gata1_data = np.array(adata.obs['GATA2/GATA1']).reshape(-1, 1)  # Reshape to 2D array
        gata2_gata1_ann = ad.AnnData(X=gata2_gata1_data)
        print(gata2_gata1_ann)

        # Update column names for the new AnnData object
        gata2_gata1_ann.obs_names = adata.obs_names
        gata2_gata1_ann.var_names = ["GATA2/GATA1"]
        print(gata2_gata1_ann)

        # Concatenate along the columns (axis=1)
        res = ad.concat([adata, gata2_gata1_ann], axis=1, merge="unique")
        print(res)
        adata = res

    # select ery cells and scale pseudotime
    adata = adata[adata.obs[CELL_CLUSTER].isin(celltypes)]
    case = ["HD", "ED BMF", "ED MDS", "ED BM normal","SDS BMF", "SDS MDS", "ED MDS/AML"]
    adata = adata[adata.obs[selection_col].isin(case)]
    print(f"Normalize 'PT_avg_combined'")
    adata.obs['PT_avg_combined'] = ensure_unit_interval(adata.obs['PT_avg_combined'])
    
    adata_hd = adata[adata.obs[selection_col].isin(["HD"])]
    adata_hd.obs['time'] = adata_hd.obs[f'PT_avg_combined']

    adata_ed_bmf = adata[adata.obs[selection_col].isin(["ED BMF"])]
    adata_ed_bmf.obs['time'] = adata_ed_bmf.obs[f'PT_avg_combined']

    adata_ed_bmn = adata[adata.obs[selection_col].isin(["ED BM normal"])]
    adata_ed_bmn.obs['time'] = adata_ed_bmn.obs[f'PT_avg_combined']

    adata_ed_mds_aml = adata[adata.obs[selection_col].isin(["ED MDS/AML"])]
    adata_ed_mds_aml.obs['time'] = adata_ed_mds_aml.obs[f'PT_avg_combined']

    adata_sds_mds = adata[adata.obs[selection_col].isin(["SDS MDS"])]
    adata_sds_mds.obs['time'] = adata_sds_mds.obs[f'PT_avg_combined']

    adata_sds_bmf = adata[adata.obs[selection_col].isin(["SDS BMF"])]
    adata_sds_bmf.obs['time'] = adata_sds_bmf.obs[f'PT_avg_combined']

    if EQ:
        print("Equalize pseudotime")
        path_base = path_base + 'eq_'
        time_min = min(max(adata_hd.obs['time']), max(adata_ed_bmf.obs['time']), max(adata_ed_bmn.obs['time']), max(adata_ed_mds_aml.obs['time']),
                    max(adata_sds_mds.obs['time']), max(adata_sds_bmf.obs['time']))
        print(f"Minim maximum time {time_min}")

        for adata in [adata_hd, adata_ed_bmf, adata_ed_bmn, adata_ed_mds_aml, adata_sds_mds, adata_sds_bmf]:
            adata = adata[adata.obs['time'] <= time_min]
        
    opt_splits = []
    for adata in [adata_hd, adata_ed_bmf, adata_ed_bmn, adata_ed_mds_aml, adata_sds_mds, adata_sds_bmf]:
        # determine optimal splitting point
        x = np.asarray(adata.obs.time)
        optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
        optb.fit(x, x)
        opt_splits.append(len(optb.splits))

    print(f"The maximum number of splits: {max(opt_splits)}")
    n_bins = max(opt_splits)

    aligner_hd_ed_bmf = Main.RefQueryAligner(adata_hd, adata_ed_bmf, marker_genes, n_bins)
    aligner_hd_ed_bmn = Main.RefQueryAligner(adata_hd, adata_ed_bmn, marker_genes, n_bins)
    aligner_hd_ed_mds_aml = Main.RefQueryAligner(adata_hd, adata_ed_mds_aml, marker_genes, n_bins)
    aligner_hd_sds_bmf = Main.RefQueryAligner(adata_hd, adata_sds_bmf, marker_genes, n_bins)
    aligner_hd_sds_mds = Main.RefQueryAligner(adata_hd, adata_sds_mds, marker_genes, n_bins)

    for aligner in [aligner_hd_ed_bmf, aligner_hd_ed_bmn, aligner_hd_ed_mds_aml, aligner_hd_sds_bmf, aligner_hd_sds_mds]:
        aligner.WEIGHT_BY_CELL_DENSITY = True
        aligner.WINDOW_SIZE=0.1
        aligner.align_all_pairs()

    # create legend
    legend_elements = [Line2D([0], [0], color='#1B9E77', label="HC",lw=3),
                    Line2D([0], [0], color='#D95F02', label="ED BMF",lw=3),
                    Line2D([0], [0], color='#D95F02', label="ED BM normal", dashes=(2, 2),lw=3),
                    Line2D([0], [0], color='#D95F02', label="ED MDS/AML", dashes=(3,2,1,2),lw=3),
                    Line2D([0], [0], color='#7570B3', label="SDS BMF",lw=3),
                    Line2D([0], [0], color='#7570B3', label="SDS MDS", dashes=(2, 2),lw=3),
                    ]
    
    df = pd.DataFrame(marker_genes).transpose()

    # prep average gene expression needed for normalizing
    hd_base= adata_hd
    hd_base = hd_base[:, marker_genes]
    hd_avg_gene_exp =  pd.DataFrame(hd_base.X.todense()) 
    hd_avg_gene_exp.columns = hd_base.var_names
    hd_avg_gene_exp = hd_avg_gene_exp.set_index(hd_base.obs_names)
    hd_avg_gene_exp =hd_avg_gene_exp.mean(axis=0)
    print(hd_avg_gene_exp)
    print(f"Mean genewise expression of HD: {hd_avg_gene_exp}")

        
    for gene in marker_genes:

        fig, ax = plt.subplots()
        al_obj_hd_ed_bmf = aligner_hd_ed_bmf.results_map[gene]
        al_obj_hd_ed_bmn = aligner_hd_ed_bmn.results_map[gene]
        al_obj_hd_ed_mds_aml = aligner_hd_ed_mds_aml.results_map[gene]
        al_obj_hd_sds_bmf = aligner_hd_sds_bmf.results_map[gene]
        al_obj_hd_sds_mds = aligner_hd_sds_mds.results_map[gene] 
        
        # scale gene expression by dividing by average hd expression 
        avg_hd = hd_avg_gene_exp[gene]

        hd_mean_trend = al_obj_hd_ed_bmf.S.mean_trend /avg_hd
        ed_bmf_mean_trend = al_obj_hd_ed_bmf.T.mean_trend/avg_hd
        ed_bmn_mean_trend = al_obj_hd_ed_bmn.T.mean_trend/avg_hd
        ed_mds_aml_mean_trend = al_obj_hd_ed_mds_aml.T.mean_trend/avg_hd
        sds_bmf_mean_trend = al_obj_hd_sds_bmf.T.mean_trend/avg_hd
        sds_mds_mean_trend = al_obj_hd_sds_mds.T.mean_trend/avg_hd

        # mean at the databins, not the interpolated mean
        ax = sb.lineplot(x= al_obj_hd_ed_bmf.S.time_points, y=hd_mean_trend, color='#1B9E77', linewidth=4) # HD
        ax = sb.lineplot(x= al_obj_hd_ed_bmf.T.time_points, y=ed_bmf_mean_trend, color='#D95F02', linewidth=4) # ED BMF
        ax = sb.lineplot(x= al_obj_hd_ed_bmn.T.time_points, y=ed_bmn_mean_trend, color='#D95F02', linewidth=4, linestyle='--') # ED BM normal
        ax = sb.lineplot(x= al_obj_hd_ed_mds_aml.T.time_points, y=ed_mds_aml_mean_trend, color='#D95F02', linewidth=4, linestyle='-.') # ED MDS / AML
        ax = sb.lineplot(x= al_obj_hd_sds_bmf.T.time_points, y=sds_bmf_mean_trend, color='#7570B3', linewidth=4) # SDS BMF
        ax = sb.lineplot(x= al_obj_hd_sds_mds.T.time_points, y=sds_mds_mean_trend, color='#7570B3', linewidth=4, linestyle='--') # SDS MDS
            
        plt.title(gene, fontsize = 20)
        plt.xlabel('Pseudotime', fontsize = 16)
        if gene == 'GATA2/GATA1':
            plt.ylabel('GATA2/GATA1', fontsize = 16)
        else:
            plt.ylabel('Gene expression (log-normalized)', fontsize = 16)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tick_params(axis='x', labelsize=14)
        plt.tick_params(axis='y', labelsize=14)

        ax.legend(handles=legend_elements, fontsize = 16,loc='upper center', bbox_to_anchor=(0.5, -0.15),  ncol=6)


        sub_path_fig = f"{path_base}/_hd_edbmf_edmds_edmdsaml_sdsbmf_sdsmds_scaled/"  #requires at the moment that the new directory name is already in the pathfig
        Path(sub_path_fig).mkdir(parents=True, exist_ok=True)

        plt.savefig(sub_path_fig + gene.replace("/", "-") + '.png', dpi=300, bbox_inches='tight')

        # add standard deviation 
        # scale standard deviation by dividing through average hd expression
        z = 1.96
        stde_hd = (al_obj_hd_ed_bmf.S.std_trend/avg_hd) / (np.sqrt(np.asarray([data_bin.shape[0] for data_bin in al_obj_hd_ed_bmf.S.data_bins])))
        stde_ed_bmf = (al_obj_hd_ed_bmf.T.std_trend/avg_hd) / (np.sqrt(np.asarray([data_bin.shape[0] for data_bin in al_obj_hd_ed_bmf.T.data_bins])))
        stde_ed_bmn = (al_obj_hd_ed_bmn.T.std_trend/avg_hd) / (np.sqrt(np.asarray([data_bin.shape[0] for data_bin in al_obj_hd_ed_bmn.T.data_bins])))
        stde_ed_mds_aml = (al_obj_hd_ed_mds_aml.T.std_trend/avg_hd) / (np.sqrt(np.asarray([data_bin.shape[0] for data_bin in al_obj_hd_ed_mds_aml.T.data_bins])))
        stde_sds_bmf = (al_obj_hd_sds_bmf.T.std_trend/avg_hd) / (np.sqrt(np.asarray([data_bin.shape[0] for data_bin in al_obj_hd_sds_bmf.T.data_bins])))
        stde_sds_mds = (al_obj_hd_sds_mds.T.std_trend/avg_hd) / (np.sqrt(np.asarray([data_bin.shape[0] for data_bin in al_obj_hd_sds_mds.T.data_bins])))
        

        ax.fill_between(al_obj_hd_ed_bmf.S.time_points, hd_mean_trend - z*stde_hd,
                             hd_mean_trend + z*stde_hd, color='#1B9E77', alpha=0.2) # HD
        ax.fill_between(al_obj_hd_ed_bmf.T.time_points, ed_bmf_mean_trend - z*stde_ed_bmf,
                             ed_bmf_mean_trend + z*stde_ed_bmf, color='#D95F02', alpha=0.2) # ED BMF
        ax.fill_between(al_obj_hd_ed_bmn.T.time_points, ed_bmn_mean_trend - z*stde_ed_bmn,
                             ed_bmn_mean_trend + z*stde_ed_bmn, color='#D95F02', alpha=0.2)# ED BM normal
        ax.fill_between(al_obj_hd_ed_mds_aml.T.time_points, ed_mds_aml_mean_trend - z*stde_ed_mds_aml,
                             ed_mds_aml_mean_trend + z*stde_ed_mds_aml, color='#D95F02', alpha=0.2) # ED MDS / AML
        ax.fill_between(al_obj_hd_sds_bmf.T.time_points, sds_bmf_mean_trend - z*stde_sds_bmf,
                             sds_bmf_mean_trend + z*stde_sds_bmf, color='#7570B3', alpha=0.2) # SDS BMF
        ax.fill_between(al_obj_hd_sds_mds.T.time_points, sds_mds_mean_trend - z*stde_sds_mds,
                             sds_mds_mean_trend + z*stde_sds_mds, color='#7570B3', alpha=0.2) # SDS MDS
        
        plt.title(gene, fontsize = 20)
        plt.xlabel('Pseudotime', fontsize = 16)
        if gene == 'GATA2/GATA1':
            plt.ylabel('GATA2/GATA1', fontsize = 16)
        else:
            plt.ylabel('Gene expression (log-normalized)', fontsize = 16)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tick_params(axis='x', labelsize=14)
        plt.tick_params(axis='y', labelsize=14)

        ax.legend(handles=legend_elements, fontsize = 16,loc='upper center', bbox_to_anchor=(0.5, -0.15),  ncol=6)

        sub_path_fig = f"{path_base}/_hd_edbmf_edmds_edmdsaml_sdsbmf_sdsmds_scaled/" 
        Path(sub_path_fig).mkdir(parents=True, exist_ok=True)

        plt.savefig(sub_path_fig + gene.replace("/", "-") + 'with_std.png', dpi=300, bbox_inches='tight')


if __name__ == "__main__": 

    # Pseutotime calculation 
    path = ""
    file_name = 'adata.h5ad'
    path_fig = "g2g/"

    store_pseudotime(path_fig, file_name, HD_CASES, HD_SELECTION_COL, path, SAVE=True)
    store_pseudotime(path_fig, file_name, ED_CASES, ED_SELECTION_COL, path, SAVE=True)
    store_pseudotime(path_fig, file_name, SDS_CASES, SDS_SELECTION_COL, path, SAVE=True)
    
    # Pseudotime bootstrapping
    file_name = 'adata_pt.h5ad'
    result_path = 'adata_bootstrapped_pt.h5ad'
    for case,col in zip(HD_CASES, HD_SELECTION_COL):
        pseudotime_bootstrapping(path_fig + 'bootstrapping_pt/', case = case, selection_col = col, 
                                 file_name= file_name, pt_file_name="g2g/pseudotimes_bootstrapped.csv")
    save_average_PT(HD_CASES, HD_SELECTION_COL, file_name, pt_file_name="g2g/pseudotimes_bootstrapped.csv", result_path=result_path)

    for case,col in zip(ED_CASES, ED_SELECTION_COL):
        pseudotime_bootstrapping(path_fig + 'bootstrapping_pt/', case = case, selection_col = col, 
                                 file_name= file_name, pt_file_name="g2g/pseudotimes_bootstrapped_ed.csv")
    save_average_PT(ED_CASES, ED_SELECTION_COL, file_name, pt_file_name="g2g/pseudotimes_bootstrapped_ed.csv", result_path=result_path)
    
    for case,col in zip(SDS_CASES, SDS_SELECTION_COL):
        pseudotime_bootstrapping(path_fig + 'bootstrapping_pt/', case = case, selection_col = col, 
                                 file_name= file_name, pt_file_name="g2g/pseudotimes_bootstrapped_sds.csv")
    save_average_PT(SDS_CASES, SDS_SELECTION_COL, file_name, pt_file_name="g2g/pseudotimes_bootstrapped_sds.csv", result_path=result_path)

    # pseudotime boostrapping combined samples
    case = ["HD", "ED BMF", "ED MDS", "ED BM normal","SDS BMF", "SDS MDS", "ED MDS/AML"]
    pseudotime_bootstrapping(path_fig + 'bootstrapping_pt/', case = case, selection_col = col, 
                             file_name= file_name, pt_file_name="g2g/pseudotimes_bootstrapped_combined_cases.csv")
    save_average_PT([case], [['condition_diagnosis']], file_name, pt_file_name="g2g/pseudotimes_bootstrapped_combined_cases.csv",
                    result_path="avg_PT_combined.h5ad")


    file_name = 'adata_bootstrapped_pt.h5ad'
    path_base = "/adjusted_by_baseline/"

    # PT alignment comparisons
    # example: HD vs ED 
    pt_comparisons(path_base, file_name, HD_REF_ED, ED_QUERY_HD, HD_ED_COL, CT_ERY, ERY_MARKER, "erythroid")
    pt_comparisons(path_base, file_name, HD_REF_ED, ED_QUERY_HD, HD_ED_COL, CT_B, B_MARKER, "b_cells")

    
    # PT alignment bootstrapping
    # example
    bootstrapping_pt_comp(path_base, file_name, HD_REF_ED, ED_QUERY_HD, HD_ED_COL, CT_ERY, ERY_MARKER)

    # PT gene alignment for three groups
    marker_genes = ['ERCC6L2', 'TP53', 'TFRC','CD36','GATA1', 'EPOR','KIT', 'ITGA4','GATA2','CD34']
    trio_comparison(path_base, file_name, celltypes_to_keep=CT_ERY, marker_genes=marker_genes)

    # PT expression plots, plots for Figure 3D
    marker_genes = ['ERCC6L2', 'TP53', 'TFRC','CD36','GATA1', 'EPOR','KIT', 'ITGA4','GATA2','CD34', 'SBDS', 'MYC', 'MYB']
    sixwise_comp(path_base, file_name, EQ=False, marker_genes=marker_genes, celltypes=CT_ERY)

