import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns
import warnings
import os
import numpy as np
import pandas as pd
from scipy import stats
import anndata as ad
import optparse, sys


cts_new = ['HSCs & MPPs', 'LMP', 'EMP',
                   'Early eryth prog', 'Late eryth prog',
                   'MkP',
                   'Granulopoietic cells',
                   'Monocytopoietic cells',
                    'DC cells',
                    'B cells',
                    'Plasma cells',
                    'T cells',
                    'NK cells',
                    'Neutrophils',
                    'Mesenchymal cells']

ct_key='celltype_short_with_neutrophils'
ct_order = ['HSCs & MPPs', 'LMP', 'EMP',
            'Early eryth prog', 'Late eryth prog', 
            'MkP', 'MkP2',
            'Early promyelo', 'Late promyelo', 'Myelo', 'EoBaso',
            'Class mono', 'Non-class mono',
            'pDC', 'cDC1', 'cDC2',
            'Small pre-B', 'Cycl pro-B & pre-B', 'Non-cycl pro-B & pre-B',
            'Immature B', 'Mature naï B', 'Non-switch mem B',
            'Plasma cells',
            'CD4+ mem T', 'CD4+ cytotoxic T', 'CD8+ naï T', 'CD8+ cent mem T',
            'NK prog', 'NK T', 'CD56bright CD16- NK',
            'Neutrophils',
            'Mesenchymal cells 1']

ct_level_ct = 'celltype_short_with_neutrophils'
cts_selected_ct = ['HSCs & MPPs', 'LMP', 'EMP',
                   'Early eryth prog', 'Late eryth prog']
ct_level_ct_coarse = 'celltype_coarse_with_neutrophils'
cts_selected_ct_coarse = ['MkP',
                          'Granulopoietic cells', 
                          'Monocytopoietic cells',
                          'DC cells',
                          'B cells',
                          'Plasma cells',
                          'T cells',
                          'NK cells',
                          'Neutrophils',
                          'Mesenchymal cells']
ct_order_selected = cts_selected_ct + cts_selected_ct_coarse


reference_cat=[
            'B cells (non-malignant)' ,
             'DC cells (non-malignant)',
            'EMP (non-malignant)',
             'Early eryth prog (non-malignant)',
             'Granulopoietic cells (non-malignant)',
             'HSCs & MPPs (non-malignant)',
             'LMP (non-malignant)',
             'Late eryth prog (non-malignant)',
             'Mesenchymal cells (non-malignant)',
             'MkP (non-malignant)',
             'Monocytopoietic cells (non-malignant)',
             'NK cells (non-malignant)',
             'Plasma cells (non-malignant)',
             'T cells (non-malignant)'
        ]



def read_data(path: str):
    """Read anndata object from the path.
    
    Args:
        path (str): filepath

    Returns:
        adata: AnnData object

    """

    adata = sc.read_h5ad(path)
    return adata


def run_infercnv(adata):
    """Run inferCNV and provide normal cell types

    Args:
        adata: AnnData object
    """

    cnv.tl.infercnv(
        adata,
        window_size=250
    )

def cluster_by_cnvprofile(adata):
    """Cluster the cells by CNV profiles with leiden clustering.
    
    Args:
        adata: AnnData object
    """

    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    cnv.tl.leiden(adata, resolution=0.35)


def umap_and_score(adata):
    """Create an umap and score the cnvs.
    """

    cnv.tl.umap(adata)
    cnv.tl.cnv_score(adata)


def plot_chromosome_heatmap(adata, groupby="cell_type", dendrogram=False):
    """Plot heatmap by chromosome (columns) and groups on rows.

    Args:
        adata: AnnData object
        groupby: grouping for rows in heatmap, eg. cell_type or cnv_cluster
        dendrogram: show dendrogram, default False
    """

    cnv.pl.chromosome_heatmap(adata, groupby=groupby, dendrogram=dendrogram)

def merge_cnv_files(all_cnvs, file):
    """Merge the full result file and sample specific file. This will be used to add CNV information to the full adata object.

    Args:
        all_cnvs, df to which all the sample specific cnv information is concatenated.
        file, df from sample adata given in the following format when calling the method: adata.obs[["cnv_leiden", "cnv_score"]]

    Returns:
        df, concatenated df with columns "cnv_leiden", "cnv_score" and cell barcodes as indexes
    """
    
    return pd.concat([all_cnvs, file])

# add gene position annotations from biomart
def annot_biomart(adata):
    """Annotate adata object with biomart gene position.
    """
    
    biomartannot = pd.read_csv("biomart_gene_post.txt", sep=" ")
    biomartannot["chromosome_name"] = "chr" + biomartannot["chromosome_name"]
    biomartannot.index = biomartannot["ensembl_gene_id"]

    #index of adata and file has to be the same
    adata.var["gene_name"] =adata.var_names
    adata.var_names = adata.var["gene_ids"]

    adata.var["chromosome"] = biomartannot["chromosome_name"]
 
    adata.var["start"] = biomartannot["start_position"]
    adata.var["end"] = biomartannot["end_position"]
    adata.var_names = adata.var["gene_name"]

def join_coarce_cell_types(adata):
    """Creates new cell type column, that groups together T-cell subtypes and same with B-cells etc.

    Returns:
	ct_key_new, the name for the new cell type column
    """

    ct_key_short = 'celltype_short_with_neutrophils'
    ct_key_coarse = 'celltype_coarse_with_neutrophils'
    ct_key_new = 'celltype_short_for_ery_coarse_else_with_neutrophils'

    observations = adata.obs[[ct_key_short, ct_key_coarse]]
    observations[ct_key_new] = observations[ct_key_coarse]

    eep = 'Early eryth prog'
    lep = 'Late eryth prog'
    indices_eep = adata[adata.obs[ct_key_short] == eep].obs.index
    indices_lep = adata[adata.obs[ct_key_short] == lep].obs.index

    observations[ct_key_new] = observations[ct_key_new].astype('str')
    observations.loc[indices_eep, ct_key_new] = eep
    observations.loc[indices_lep, ct_key_new] = lep
    observations[ct_key_new] = observations[ct_key_new].astype('category')

    adata.obs[ct_key_new] = observations[ct_key_new]

    return ct_key_new

def reorder_cell_types(adata, cts_new, ct_key_new):

    adata.obs[ct_key_new] = adata.obs[ct_key_new].cat.reorder_categories(cts_new)

def plot_cnv_umap(adata,title):
    sc.settings.set_figure_params(fontsize=20)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
    ax4.axis("off")
    cnv.pl.umap(
        adata,
        color="cnv_leiden",
        legend_loc="on data",
        legend_fontoutline=2,
        ax=ax1,
        show=False,
    )
    cnv.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
    cnv.pl.umap(adata, color=ct_key, ax=ax3)
    fig.savefig(title, dpi=300,bbox_inches="tight") ### removed path os.path.join(path, title)


def plot_laura_umap(adata,title):
    sc.set_figure_params(fontsize=20)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"wspace": 0.5})
    ax4.axis("off")
    sc.pl.umap(adata, color="cnv_leiden", ax=ax1, show=False)
    sc.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
    sc.pl.umap(adata, color=ct_key, ax=ax3)
    
    fig.savefig(title, dpi=300,bbox_inches="tight")

def plot_single_umap(adata, color, filename):

    with rc_context({"figure.figsize": (8,8)}): 
        sc.pl.umap(adata, color=color, size=40, save=filename, frameon=False, title="CNV score")


def plot_single_sample_umap(adata, sample, color, filename):

    with rc_context({"figure.figsize": (8,8)}): 
        sc.pl.umap(adata, groups=[sample], color=color, size=40, save=filename, frameon=False, title=f"CNV score") ##, {sample}
        

def plot_multiple_single_sample_umap(adata, samples, color, filename):

    for sample in samples:
        test = adata.copy()
        test.obs["cnv_score"] = np.nan
        test.obs["cnv_score"] = adata[adata.obs["sample"]==sample].obs["cnv_score"]
        plot_single_sample_umap(test, sample=sample, color="cnv_score", filename=f"{sample}_{filename}.png")


def violinplot(adata,filename, groupby='celltype_short_for_ery_coarse_else_with_neutrophils', ylabel="cnv_score", order=cts_new, stripplot=False):

    sc.set_figure_params(fontsize=40)
    
    with rc_context({"figure.figsize": (15,15)}): 
        ax = sc.pl.violin(adata,
                ["cnv_score"],
                groupby=ct_key_new,
                stripplot=False,
                    rotation=90,
                         show=False,
                         order=order)
        plt.title(f"{sample}")
        plt.savefig(filename, dpi=300, bbox_inches="tight")


def proportions_for_barplot(adata):
    sizes = adata.obs.groupby([ct_key_new, "cnv_leiden"]).size()
    props = pd.DataFrame(sizes.groupby(ct_key_new).apply(lambda x: 100 * x/ x.sum()).reset_index(level=0, drop=True).reset_index())
    props = props.pivot(columns=ct_key_new, index="cnv_leiden").T
    props.index = props.index.droplevel(0)

    return props

def proportions_gt_for_barplot(adata):
    sizes = adata.obs.groupby(["cnv_leiden", "tp53_genotype_wtmut_noTP53clone_wt"]).size()
    props = pd.DataFrame(sizes.groupby("cnv_leiden").apply(lambda x: 100 * x/ x.sum()).reset_index(level=0, drop=True).reset_index())
    props = props.pivot(columns="cnv_leiden", index="tp53_genotype_wtmut_noTP53clone_wt").T
    props.index = props.index.droplevel(0)

    return props

def plot_cluster_proportions(cluster_props, filename,title,
                             cluster_palette=None,
                             xlabel_rotation=0,
                            figsize=(8,8)
                            ): 
    sc.set_figure_params(fontsize=30)
    fig, ax = plt.subplots(dpi=300, figsize=figsize)
    fig.patch.set_facecolor("white")
    
    cmap = None
    if cluster_palette is not None:
        cmap = sns.palettes.blend_palette(
            cluster_palette, 
            n_colors=len(cluster_palette), 
            as_cmap=True)
   
    cluster_props.plot(
        kind="bar", 
        stacked=True, 
        ax=ax, 
        legend=None, 
        colormap=cmap
    )

    #reorder legend labels
    handles, labels = ax.get_legend_handles_labels()
    labels = list(reversed(labels))
    handles = list(reversed(handles))
    ax.legend(handles, labels, bbox_to_anchor=(1.01, 1), frameon=False, title=title)
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=xlabel_rotation)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Proportion")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    
    return fig

def get_ct_gene(adata, sample, gene, path, ct_key_new):

    gene_adata = adata[:,adata.var_names ==gene]
    results = pd.DataFrame(columns=["slope", "intercept", "r", "p", "std_err"])

    for ct in set(adata.obs[ct_key_new]):

        try:
            slope, intercept, r, p, std_err = gene_expression_cnvScore(gene_adata[gene_adata.obs[ct_key_new]==ct], gene, ct, sample)
            results.loc[ct, [ "slope", "intercept", "r", "p", "std_err"]] = [slope, intercept, r, p, std_err ]
        except:
            results.loc[ct, [ "slope", "intercept", "r", "p", "std_err"]] = [np.nan, np.nan, np.nan, np.nan, np.nan]

    results.to_csv(os.path.join(path,f"{sample}_{gene}_expr_cnvscore_lm.txt"), sep="\t")


def gene_expression_cnvScore(adata, gene, cell, sample):
    """Linear regression gene_expression ~ cnv_score.
    """
    gene_expr = adata.X.toarray()
    cnv_score = adata.obs["cnv_score"]
    try:
        slope, intercept, r, p, std_err = stats.linregress(gene_expr[:,0], cnv_score)
    except:
    	return np.nan
    plt.scatter(cnv_score, gene_expr[:,0])
    plt.ylabel(f"{gene} expression")
    plt.xlabel("cnv score")
    plt.title(cell)
    plt.savefig(os.path.join(path,f"{sample}_{cell}_{gene}expression_cnvScore_scatter.png"),dpi=300)
    plt.close()

    return (slope, intercept, r, p,std_err)

def merge_cnv_files(all_cnvs, file):
    """Merge the full result file and sample specific file. This will be used to add CNV information to the full adata object.

    Args:
	all_cnvs, df to which all the sample specific cnv information is concatenated.
        file, df from sample adata given in the following format when calling the method: adata.obs[["cnv_leiden", "cnv_score"]]

    Returns:
	df, concatenated df with columns "cnv_leiden", "cnv_score" and cell barcodes as indexes
    """

    return pd.concat([all_cnvs, file])


def write_adata(adata, filename):

    adata.write_h5ad(filename)

    print(f"Saved file {filename}.")


def set_control_celltypes(adata):
    temp = pd.DataFrame(adata.obs[["sample",ct_key_new]])

    ct_annotation = {'B cells': 'B cells (non-malignant)' ,
     'DC cells': 'DC cells (non-malignant)',
     'EMP': 'EMP (non-malignant)',
     'Early eryth prog': 'Early eryth prog (non-malignant)',
     'Granulopoietic cells': 'Granulopoietic cells (non-malignant)',
     'HSCs & MPPs': 'HSCs & MPPs (non-malignant)',
     'LMP': 'LMP (non-malignant)',
     'Late eryth prog': 'Late eryth prog (non-malignant)',
     'Mesenchymal cells': 'Mesenchymal cells (non-malignant)',
     'MkP': 'MkP (non-malignant)',
     'Monocytopoietic cells': 'Monocytopoietic cells (non-malignant)',
     'NK cells': 'NK cells (non-malignant)',
     'Neutrophils': 'Neutrophils (non-malignant)',
     'Plasma cells': 'Plasma cells (non-malignant)',
     'T cells': 'T cells (non-malignant)'}
    
    temp["ct_cnv"] = pd.DataFrame(temp.loc[temp["sample"].isin(healthy_samples), ct_key_new].map(ct_annotation))
    temp["test"]= (temp["ct_cnv"]== np.nan).astype(str)
    temp["ct_cnv"] = temp["ct_cnv"].astype(str, errors="ignore")
    temp.loc[temp["ct_cnv"]=="nan", "ct_cnv"] = temp.loc[temp["ct_cnv"]=="nan", ct_key_new]
    
    adata.obs["ct_cnv"] = temp["ct_cnv"]


def visuals(adata, adata_orig, sample):

    cnv.pl.chromosome_heatmap(adata, groupby=ct_key_new, figsize=(15,23), save=f"{sample}_chromosome_heatmapt_ct.png")
    cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", save=f"{sample}_chromosome_heatmapt_leiden.png")

    # reorder doesn't work if the file doesn't have all the cell types listed!
    reorder_cell_types(adata, ct_order_selected, ct_key_new)
    plot_cnv_umap(adata,f"{sample}_vsControl_perct_cnv_UMAPs_leiden_cnvscore.png")
    plot_laura_umap(adata, f"{sample}_vsControl_perct_scanpy_UMAPs_leiden_cnvscore.png")
    violinplot(adata, groupby=ct_key_new, ylabel="cnv_score", order=ct_order_selected, stripplot=False,filename=f"{sample}_vsControl_perct_cnvscore_ctcoarse_violinplot.png")

    adata_orig.obs["cnv_score"] = adata[adata.obs["sample"]==sample].obs["cnv_score"]
    plot_multiple_single_sample_umap(adata_orig, [sample], "cnv_score", "_cnvscore_umap")
    try:
        reorder_cell_types(adata, cts_new, ct_key_new)
    except:
        print("Couldn't reorder cells")
    plot_cluster_proportions(proportions_for_barplot(adata),filename=f"{sample}_proportion_ct.png", title='CNV\ncluster',xlabel_rotation=90, figsize=(10,10) )
    try:
        plot_cluster_proportions(proportions_for_barplot(adata).loc[["HSCs & MPPs", "EMP", "Early eryth prog", "Late eryth prog"],:],filename=f"{sample}_proportion_ct_HSCT_EEP_LEP.png", title='CNV\ncluster',xlabel_rotation=90, figsize=(8,8))
    except:
        plot_cluster_proportions(proportions_for_barplot(adata).loc[["EMP", "Early eryth prog", "Late eryth prog"],:],filename=f"{sample}_proportion_ct_HSCT_EEP_LEP.png", title='CNV\ncluster',xlabel_rotation=90, figsize=(8,8))
    try:    
        plot_cluster_proportions(proportions_gt_for_barplot(adata), filename=f"{sample}_proportion_TP53.png", title='TP53\ngenotype',xlabel_rotation=90, figsize=(10,10) )
    except:
        print(f"Couldn't plot TP53 proportions")


def infercnv_per_ct(subset, adata_orig, sample, reference_samples):
    ads = []
    
    i=0
    
    for ct in set(subset.obs[ct_key_new]):
        print(ct)
        print(ct_annotation[ct])

        tmp = subset[subset.obs[ct_key_new].isin([ct])]
        try:
            if ct == "Neutrophils":
                 new=cnv.tl.infercnv(
                    tmp,
                    window_size=250,
                     reference_samples=reference_samples,
                )
            else:
    
                new=cnv.tl.infercnv(
                    tmp,
                    reference_key="ct_ref", #all controls as reference!
                    reference_cat=[*ct_annotation[ct]], 
                    reference_samples=reference_samples,
                    window_size=250
                )
            
            ads.append(new)
        except:
            print(f'Sample has no {ct} cells.')
    
    concat = ad.concat(ads,merge="same", uns_merge="same")

    cluster_by_cnvprofile(concat)
    umap_and_score(concat)

    visuals(concat, adata_orig, sample)

    return concat


def make_control(n: int, healthy_samples: list):
    """Get a random set of healthy samples to be used as a control.
        Args:
            n, number of samples in control
            healthy_samples, list of names of the healthy samples

        Returns:
            control_sample, list of names of the samples in control
    """

    import random
    random.seed(42)


    control_sample = random.sample(healthy_samples, i)
    return control_sample

def create_ctref(adata):
    """Make coarse ct_key and add ct_ref annotation for the control cells."""

    ct_key_new = join_coarce_cell_types(adata)
    set_control_celltypes(adata)
    #takes only control samples
    temp = pd.DataFrame(adata[adata.obs["ct_cnv"].isin(reference_cat)].obs[["sample","ct_cnv",ct_key_new]])
    temp2 = temp.loc[~temp["sample"].isin(control_sample)]
    temp2["ct_ref"] = temp2[ct_key_new]
    temp = temp.loc[temp["sample"].isin(control_sample)]
    temp["ct_ref"] = temp["sample"].astype(str) + " " + temp["ct_cnv"].astype(str)
    temp["ct_ref"] = temp["ct_ref"].astype(str, errors="ignore")

    adata.obs["ct_ref"] = temp["ct_ref"]

    temp = temp.groupby(ct_key_new).agg({"ct_ref": lambda x: list(x.unique())}).reset_index()
    ct_annotation = dict(zip(temp[ct_key_new],temp["ct_ref"]))
    ct_annotation['Neutrophils'] = list(reference_cat)

    return ct_key_new, ct_annotation


def main(args):


    adata = read_data(os.path.join(args.folder, args.file))
    control_sample = make_control(args.n, args.healthy_samples)
    ct_key_new, ct_annotation = create_ctref(adata)


    annot_biomart(adata)
    adata = clean_adata(adata)

    for sample in args.samples:
        samples = control_sample + [sample]
        sample_adata = adata[atada.obs["sample"].isin(samples)]

        adata_infercnv = infercnv_per_ct(sample_adata, adata, sample, control_sample)


if __name__ == "__main__":

    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", dest="input", action="store", help="The input anndata file", metavar="data/adata.h5ad")
    parser.add_option("-f", "--folder", dest="folder", action="store", help="The path to folder where the input file is")
    parser.add_option("-n", "--nsamples", dest="nsamples", action="store", help="Number of healthy samples used in a control set")
    parser.add_option("-hs", "--hsamples", dest="hsamples", action="store", help="List of names of healthy samples")

    (options, args) = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(options)
