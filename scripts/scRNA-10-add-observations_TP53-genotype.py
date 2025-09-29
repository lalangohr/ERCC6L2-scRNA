"""
Part 10 of scRNA-seq data analysis. 

Add TP53 genotype (WT/MUT) to adata object.
"""

import scanpy as sc
import pandas as pd
import csv
import os


# Folder in which result is written
dir_results = './results'


adata_file = 'adata_filtered_seeds_scANVI_UMAP_unclear-cell-types_cell-cycle.h5ad'
adata = sc.read_h5ad(os.path.join(dir_results, adata_file))


# Classes

class Primer_stats:
    """
    Class for primer statistics for one variant
    """

    individual = ''
    variant = ''
    primer = ''

    n_cells = 0
    primer_counts_fractions = pd.DataFrame()

    def calculate_stats_for_one_variant(self):

        adata_subset = adata[adata.obs['sample'].isin([self.individual])]

        self.n_cells = adata_subset.n_obs

        primer_counts = adata_subset.obs.groupby([self.individual + '_' + self.variant + '_' + self.primer], dropna=False).apply(len)
        
        primer_fractions = primer_counts / adata_subset.n_obs

        primer_counts_fractions = pd.concat([primer_counts, primer_fractions], axis=1, join="inner")

        self.primer_counts_fractions = primer_counts_fractions.set_axis(['count','fraction'], axis=1)

        print(f"\n#cells: {self.n_cells}")
        print(f"\n{self.primer_counts_fractions}")


class Primer_stats_all():
    """
    Class for all primer statistics
    """

    list_stats = []

    def add_primer_stats(self, curr_primer_stats, individual, variant, primer):

        if len(primer) == 1:
            primer_string = primer[0]
        else:
            primer_string = primer[0] + ' & ' + primer[1]

        list_stats_curr = [individual, variant, primer, primer_string]
        list_stats_curr.append(curr_primer_stats.n_cells)
        list_stats_curr.extend(curr_primer_stats.primer_counts_fractions[['count','fraction']].sum().tolist())
        self.list_stats.append(list_stats_curr)

        print(f"list_stats_curr: {list_stats_curr}")

    def write_to_file(self):

        df_stats = pd.DataFrame(self.list_stats,
                                columns =['individual', 'variant', 'primer', 'primer_string', '#all', '#genotyped', '%genotyped'])
        df_stats = df_stats.astype({"#all": int, "#genotyped": int})
        df_stats['%genotyped'] = df_stats['%genotyped'].apply(lambda x: x*100)
        df_stats.to_csv(dir_results+'tp53_genotype_stats_all.csv', index=False)


class Two_primer_stats():
    """
    Class for primer statistics for one variant where data for two primers exits
    """

    n_all = 0
    n_identical = 0
    n_conflicting = 0
    n_conflicting_in_adata = 0


class Two_primer_stats_all():
    """
    Class for primer statistics where data for two primers exits
    """

    list_stats = []

    def add_two_primer_stats(self, current_two_primer_stats, individual, variant):

        current_list = [
            individual, 
            variant, 
            current_two_primer_stats.n_all, 
            current_two_primer_stats.n_identical, 
            current_two_primer_stats.n_conflicting, 
            current_two_primer_stats.n_conflicting_in_adata]
        self.list_stats.append(current_list)

    def write_to_file(self):

        df_stats = pd.DataFrame(self.list_stats,
                                columns =['individual', 'variant', 'n_all', 'n_identical', 'n_conflicting', 'n_conflicting_in_adata'])
        df_stats.to_csv(dir_results+'tp53_genotype_two_primer_stats_all.csv', index=False)


# Functions

def get_tp53_files(type, method):
    """
    Get files where TP53 genotype (WT/MUT) information was saved to in scAmp-seq part 8.
    """

    type = '_' + type
    method = '_' + method

    bc_filename = 'barcodes' + method + type + '.csv'

    list_individuals = [
        'E1a',
        'E1b',
        'E1c',
        # E2 has no TP53 clone
        'E3',
        'E4a', 'E4a', 'E4a',
        'E4b', 'E4b', 'E4b', 'E4b',
        'E4c', 'E4c', 'E4c', 'E4c',
        'E5', 'E5', 'E5', 'E5',
        'E6', 'E6',
        'E8a', 'E8a',
        'E8b', 'E8b', 'E8b', 'E8b',
        'E8c', 'E8c', 'E8c', 'E8c',
        'E7', 'E7',
        'E9',
        'E10',
        # SDS1 has no TP53 clone
        # SDS2 has no TP53 clone
        'SDS3',
        'SDS4', 'SDS4', 'SDS4',
        'SDS5a',
        'SDS5b'
    ]

    list_variants = [
        'c.490A>G',                                                  # E1a
        'c.490A>G',                                                  # E1b
        'c.490A>G',                                                  # E1c
                                                                     # E2 has no TP53 clone
        'c.524G>A',                                                  # E3
        'c.743G>A', 'c.818G>A', 'c.814G>A',                          # E4a
        'c.743G>A', 'c.818G>A', 'c.814G>A', 'c.41T>C',               # E4b
        'c.743G>A', 'c.818G>A', 'c.814G>A', 'c.41T>C',               # E4c
        'c.743G>A', 'c.830G>T', 'c.528C>G', 'c.843C>A',              # E5
        'c.723delC', 'c.733G>A',                                     # E6
        'c.659A>G', 'c.725G>T',                                      # E8a
        'c.659A>G', 'c.725G>T', 'c.764T>G', 'c.795_796insCCTCTTGCT', # E8b
        'c.659A>G', 'c.725G>T', 'c.764T>G', 'c.795_796insCCTCTTGCT', # E8c
        'c.659A>G', 'c.716A>G',                                      # E7
        'c.713G>A',                                                  # E9
        'c.817C>T',                                                  # E10
                                                                     # SDS1 has no TP53 clone
                                                                     # SDS2 has no TP53 clone
        'c.524G>A',                                                  # SDS3
        'c.329G>C', 'c.818G>A', 'c.833C>T',                          # SDS4
        'c.713G>C',                                                  # SDS5a
        'c.713G>C'                                                   # SDS5b
    ]
    
    list_primers = [
        ['scAmp3'],                                                                 # E1a
        ['scAmp3'],                                                                 # E1b
        ['scAmp3'],                                                                 # E1c
                                                                                    # E2 has no TP53 clone
        ['scAmp3'],                                                                 # E3
        ['scAmp4'], ['scAmp5', 'scAmp5_II'], ['scAmp5', 'scAmp5_II'],               # E4a
        ['scAmp4'], ['scAmp5', 'scAmp5_II'], ['scAmp5', 'scAmp5_II'], ['scAmp1'],   # E4b
        ['scAmp4'], ['scAmp5', 'scAmp5_II'], ['scAmp5', 'scAmp5_II'], ['scAmp1'],   # E4c
        ['scAmp4'], ['scAmp5', 'scAmp5_II'], ['scAmp3'], ['scAmp5', 'scAmp5_II'],   # E5
        ['scAmp4'], ['scAmp4'],                                                     # E6
        ['scAmp4'], ['scAmp4'],                                                     # E7a
        ['scAmp4'], ['scAmp4'], ['scAmp4', 'scAmp5'], ['scAmp5', 'scAmp5_II'],      # E7b
        ['scAmp4'], ['scAmp4'], ['scAmp4', 'scAmp5'], ['scAmp5', 'scAmp5_II'],      # E7c
        ['scAmp4'], ['scAmp4'],                                                     # E8
        ['scAmp4'],                                                                 # E9
        ['scAmp5', 'scAmp5_II'],                                                    # E10
                                                                                    # SDS1 has no TP53 clone
                                                                                    # SDS2 has no TP53 clone
        ['scAmp3'],                                                                 # SDS3
        ['scAmp2'], ['scAmp5', 'scAmp5_II'], ['scAmp5', 'scAmp5_II'],               # SDS4
        ['scAmp4'],                                                                 # SDS5a
        ['scAmp4']                                                                  # SDS5b
    ]

    list_files = [
        ['./results/E1a_c490AG/' + bc_filename], # E1a
        ['./results/E1b_c490AG/' + bc_filename], # E1b
        ['./results/E1c_c490AG/' + bc_filename], # E1c
        # E2 has no TP53 clone
        ['./results/E3_c524GA/' + bc_filename], # E3
        ['./results/E4a_c743GA/' + bc_filename], # E4a
        ['./results/E4a_c818GA/' + bc_filename, './results/E4a_c818GA-scAmp5_II/' + bc_filename], # E4a
        ['./results/E4a_c814GA/' + bc_filename, './results/E4a_c814GA-scAmp5_II/' + bc_filename], # E4a
        ['./results/E4b_c743GA/' + bc_filename], # E4b
        ['./results/E4b_c818GA/' + bc_filename, './results/E4b_c818GA-scAmp5_II/' + bc_filename], # E4b
        ['./results/E4b_c814GA/' + bc_filename, './results/E4b_c814GA-scAmp5_II/' + bc_filename], # E4b
        ['./results/E4b_c41TC/' + bc_filename],  # E4b
        ['./results/E4c_c743GA/' + bc_filename], # E4c
        ['./results/E4c_c818GA/' + bc_filename, './results/E4c_c818GA-scAmp5_II/' + bc_filename], # E4c
        ['./results/E4c_c814GA/' + bc_filename, './results/E4c_c814GA-scAmp5_II/' + bc_filename], # E4c
        ['./results/E4c_c41TC/' + bc_filename], # E4c
        ['./results/E5_c743GA/' + bc_filename], # E5
        ['./results/E5_c830GT/' + bc_filename, './results/E5_c830GT-scAmp5_II/' + bc_filename], # E5
        ['./results/E5_c528CG/' + bc_filename], # E5
        ['./results/E5_c843CA/' + bc_filename, './results/E5_c843CA-scAmp5_II/' + bc_filename], # E5
        ['./results/E6_c723delC/' + bc_filename], # E6
        ['./results/E6_c733GA/' + bc_filename],  # E6
        ['./results/E8a_c659AG/' + bc_filename], # E8a
        ['./results/E8a_c725GT/' + bc_filename], # E8a
        ['./results/E8b_c659AG/' + bc_filename], # E8b
        ['./results/E8b_c725GT/' + bc_filename], # E8b
        ['./results/E8b_c764TG-scAmp4/' + bc_filename, './results/E8b_c764TG-scAmp5/' + bc_filename], # E8b
        ['./results/E8b_c795_796insCCTCTTGCT/' + bc_filename, './results/E8b_c795_796insCCTCTTGCT-scAmp5_II/' + bc_filename], # E8b
        ['./results/E8c_c659AG/' + bc_filename],
        ['./results/E8c_c725GT/' + bc_filename],
        ['./results/E8c_c764TG-scAmp4/' + bc_filename, './results/E8c_c764TG-scAmp5/' + bc_filename],
        ['./results/E8c_c795_796insCCTCTTGCT/' + bc_filename, './results/E8c_c795_796insCCTCTTGCT-scAmp5_II/' + bc_filename],
        ['./results/E7_c659AG/' + bc_filename], # E7
        ['./results/E7_c716AG/' + bc_filename], # E7
        ['./results/E9_c713GA/' + bc_filename], # E9
        ['./results/E10_c817CT/' + bc_filename, './results/E10_c817CT-scAmp5_II/' + bc_filename], # E10
        ['./results/SDS3_c524GA/' + bc_filename], # SDS3
        ['./results/SDS4_c329GC/' + bc_filename], # SDS4
        ['./results/SDS4_c818GA/' + bc_filename, './results/SDS4_c818GA-scAmp5_II/' + bc_filename], # SDS4
        ['./results/SDS4_c833CT/' + bc_filename, './results/SDS4_c833CT-scAmp5_II/' + bc_filename], # SDS4
        ['./results/SDS5a_c713GC/' + bc_filename], # SDS5a
        ['./results/SDS5b_c713GC/' + bc_filename]  # SDS5b
    ]

    df_tp53_files = pd.DataFrame({
        'individual': list_individuals,
        'variant'   : list_variants,
        'primer'    : list_primers,
        'file'      : list_files 
    })

    return df_tp53_files


def get_bcs_from_adata():
    """
    Get barcode (bs) from adata
    """

    list_bc = adata.obs_names.tolist()
    
    des_bc = ['barcode']
    
    df_bc = pd.DataFrame(list_bc, columns=des_bc) 
    df_bc = df_bc.set_index('barcode')
    
    return df_bc


def parse_gt(file, individual, variant, suffix_capital ):
    """
    Parse TP53 genotype (gt) from file
    """
    
    df_tp53_ind_var = pd.read_csv(file)
    df_tp53_ind_var = df_tp53_ind_var.set_index('barcode')
    df_tp53_ind_var.index += '-' + individual

    df_tp53_ind_var = df_tp53_ind_var.replace({'TP53_gt_' + suffix_capital : '0'}, 'WT')
    df_tp53_ind_var = df_tp53_ind_var.replace({'TP53_gt_' + suffix_capital : '1'}, 'MUT')

    return df_tp53_ind_var


def handle_duplicates(df_tp53_ind_var_bothPrimers, individual, variant, df_bc, two_primer_stats_all):
    """
    Handle duplicates, i.e. barcodes that were genotyped in both primer files
    """

    current_two_primer_stats = Two_primer_stats()

    # Duplicate indices, i.e. barcodes which are genotyped in both primer files 
    df_duplicated = df_tp53_ind_var_bothPrimers[df_tp53_ind_var_bothPrimers.index.duplicated(keep=False)]
    current_two_primer_stats.n_all = int(len(df_duplicated.index)/2)
    
    df_duplicated = df_duplicated.reset_index(level=0)
    current_two_primer_stats.n_identical = len(df_duplicated[df_duplicated.duplicated()])
    
    # Identify duplicate indices for which genotypes differ
    df_duplicated_differ = df_duplicated[-df_duplicated.duplicated(keep=False)].sort_values(by=['barcode']).set_index('barcode')
    current_two_primer_stats.n_conflicting = int(len(df_duplicated_differ.index)/2)

    # Remove rows with duplicate indices for which genotypes differ
    df_tp53_ind_var_bothPrimers = df_tp53_ind_var_bothPrimers.drop(df_duplicated_differ.index.tolist())

    # For the remaining, remove rows with duplicate indices: keep first only
    df_tp53_ind_var_bothPrimers = df_tp53_ind_var_bothPrimers.loc[~df_tp53_ind_var_bothPrimers.index.duplicated(keep='first')]

    # Identify how many of the duplicate indices for which genotypes differ does actually appear in the scRNA-seq data (i.e. anndata object)
    set_dup_idx = df_duplicated_differ.index
    set_idx = df_bc.index
    set_int = set_idx.intersection(set_dup_idx)
    current_two_primer_stats.n_conflicting_in_adata = len(set_int)

    two_primer_stats_all.add_two_primer_stats(current_two_primer_stats, individual, variant)

    return df_tp53_ind_var_bothPrimers


def add_gt_of_one_variant_to_adata(individual, variant, primers, files, suffix, df_bc, two_primer_stats_all): 
    """
    Add TP53 genotype (gt) information of one variant to adata.
    """

    if suffix == 'wtmut':
        suffix_capital = 'WTMUT'
    elif suffix == 'wthethom':
        suffix_capital = 'WTHETHOM'

    if len(files) == 1:

        # Parse genotype info
        df_tp53_ind_var = parse_gt(files[0], individual, variant, suffix_capital)

        # Make sure to include only barcodes present in scRNA-seq data
        df_bc_tp53_ind_var = pd.concat([df_bc, df_tp53_ind_var], axis=1, join="inner").reindex(df_bc.index)

        # Add genotype as observation to the anndata object
        adata.obs[individual + '_' + variant + '_' + primers[0]] = pd.Categorical(df_bc_tp53_ind_var['TP53_gt_' + suffix_capital])


    elif len(files) == 2:

        # Parse genotype info from 1st file
        df_tp53_ind_var_first_primer = parse_gt(files[0], individual, variant, suffix_capital)

        # Parse genotype info from 2nd file
        df_tp53_ind_var_second_primer = parse_gt(files[1], individual, variant, suffix_capital)

        # Append genotype info from 1st and 2nd file
        df_tp53_ind_var_both_primers = pd.concat([df_tp53_ind_var_first_primer, df_tp53_ind_var_second_primer])

        df_tp53_ind_var_both_primers = handle_duplicates(df_tp53_ind_var_both_primers, individual, variant, df_bc, two_primer_stats_all)

        # Make sure to include only barcodes present in scRNA-seq data
        df_bc_tp53_ind_var_first_primer = pd.concat([df_bc, df_tp53_ind_var_first_primer], axis=1, join="inner").reindex(df_bc.index)
        df_bc_tp53_ind_var_second_primer = pd.concat([df_bc, df_tp53_ind_var_second_primer], axis=1, join="inner").reindex(df_bc.index)
        df_bc_tp53_ind_var_both_primers = pd.concat([df_bc, df_tp53_ind_var_both_primers], axis=1, join="inner").reindex(df_bc.index)

        # Add genotype as observation to anndata object
        adata.obs[individual + '_' + variant + '_' + primers[0]] = pd.Categorical(df_bc_tp53_ind_var_first_primer['TP53_gt_' + suffix_capital])
        adata.obs[individual + '_' + variant + '_' + primers[1]] = pd.Categorical(df_bc_tp53_ind_var_second_primer['TP53_gt_' + suffix_capital])
        adata.obs[individual + '_' + variant + '_' + primers[0] + '_' + primers[1]] = pd.Categorical(df_bc_tp53_ind_var_both_primers['TP53_gt_' + suffix_capital])


def add_gt_to_adata(df_tp53_files, suffix, df_bc, two_primer_stats_all):
    """
    Add TP53 genotype (gt) information to adata.
    """

    all_primer_stats = Primer_stats_all()

    for index, df_index in df_tp53_files.T.items():

        individual = df_index['individual']
        variant = df_index['variant']
        primers = df_index['primer']
        files = df_index['file']

        add_gt_of_one_variant_to_adata(individual, variant, primers, files, suffix, df_bc, two_primer_stats_all)

        # Statistics
        current_primer_stats = Primer_stats()
        current_primer_stats.individual = individual
        current_primer_stats.variant = variant
        if len(primers) == 1:
            current_primer_stats.primer = primers[0]
            current_primer_stats.calculate_stats_for_one_variant()
            all_primer_stats.add_primer_stats(current_primer_stats, individual, variant, primers)
        if len(primers) == 2:
            current_primer_stats.primer = primers[0]
            current_primer_stats.calculate_stats_for_one_variant()
            all_primer_stats.add_primer_stats(current_primer_stats, individual, variant, [primers[0]])

            current_primer_stats.primer = primers[1]
            current_primer_stats.calculate_stats_for_one_variant()
            all_primer_stats.add_primer_stats(current_primer_stats, individual, variant, [primers[1]])

            current_primer_stats.primer = primers[0] + '_' + primers[1]
            current_primer_stats.calculate_stats_for_one_variant()
            all_primer_stats.add_primer_stats(current_primer_stats, individual, variant, primers)

    all_primer_stats.write_to_file()


def add_overall_genotype(df_tp53_files, suffix, individuals_without_clone):
    """
    Combine genotype of different variants into overall genotype.
    """

    cell_ids_mut = pd.Index([])
    cell_ids_wt = pd.Index([])

    for individual, df_ind in df_tp53_files.groupby('individual'):

        cell_ids_ind_mut = pd.Index([])
        cell_ids_ind_wt = (adata.obs).index

        # Iterate over all variants of current individual
        for row_index, row in df_ind.iterrows():
            variant = row['variant']
            primers = row['primer']
            if len(primers) == 1:
                current_primer = primers[0]
            if len(primers) == 2:
                current_primer = primers[0] + '_' + primers[1]

            cell_ids_ind_var_mut = (
                adata.obs
                    .loc[lambda x: x['sample'] == individual]
                    .loc[lambda x: x[individual + '_' + variant + '_' + current_primer] == 'MUT'] # primer
            ).index
            cell_ids_ind_var_wt = (
                adata.obs
                    .loc[lambda x: x['sample'] == individual]
                    .loc[lambda x: x[individual + '_' + variant + '_' + current_primer] == 'WT'] # primer
            ).index

            # If at least one variant location is mutated, then barcode is defined as mutated
            cell_ids_ind_mut = cell_ids_ind_mut.union(cell_ids_ind_var_mut)

            # If all variant locations are wild-type, then barcode is defined as wild-type
            cell_ids_ind_wt = cell_ids_ind_wt.intersection(cell_ids_ind_var_wt)

        # Combine genotype info of all individuals
        cell_ids_mut = cell_ids_mut.union(cell_ids_ind_mut)
        cell_ids_wt = cell_ids_wt.union(cell_ids_ind_wt)


    # Add genotype as observation to the anndata object
    adata.obs.loc[cell_ids_mut, 'tp53_genotype_' + suffix] = 'MUT'
    adata.obs.loc[cell_ids_wt, 'tp53_genotype_' + suffix] = 'WT'


    # Mark all cells from patients without TP53 clones as WT
    cell_ids_noTP53clone = (adata[adata.obs['individual'].isin(individuals_without_clone)].obs).index

    # Add genotype as observation to the anndata object
    adata.obs.loc[cell_ids_mut, 'tp53_genotype_' + suffix + '_noTP53clone_wt'] = 'MUT'
    adata.obs.loc[cell_ids_wt, 'tp53_genotype_' + suffix + '_noTP53clone_wt'] = 'WT'
    adata.obs.loc[cell_ids_noTP53clone, 'tp53_genotype_' + suffix + '_noTP53clone_wt'] = 'WT'


def write_overall_stats(suffix):
    """
    Write statistics to file.
    """

    # Get subset of ERCC6L2 syndrome samples only to calculate some statistics
    adata_subset = adata[adata.obs['condition'].isin(['ERCC6L2 disease', 'SDS'])]

    # Calculate some statistics - over all ERCC6L2 and SDS patients
    all_counts = adata_subset.obs.groupby(['tp53_genotype_' + suffix], dropna=False).apply(len)
    all_fractions = adata_subset.obs.groupby(['tp53_genotype_' + suffix], dropna=False).apply(len) / adata_subset.n_obs
    all_counts_fractions = pd.concat([all_counts, all_fractions], axis=1, join="inner")
    all_counts_fractions = all_counts_fractions.set_axis(['count', 'fraction'], axis=1)

    all_counts_fractions.to_csv(dir_results+'tp53_genotype_' + suffix + '_overall.csv',na_rep='NaN')

    # Calculate some statistics - per patient
    ind_counts = adata_subset.obs.groupby(['sample', 'tp53_genotype_' + suffix], dropna=False).apply(len)
    ind_fractions = adata_subset.obs.groupby(['sample', 'tp53_genotype_' + suffix], dropna=False).apply(len) / adata_subset.obs.groupby(['sample'], dropna=False).apply(len)
    ind_counts_fractions = pd.concat([ind_counts, ind_fractions], axis=1, join="inner")
    ind_counts_fractions = ind_counts_fractions.set_axis(['count','fraction'], axis=1)

    ind_counts_fractions.to_csv(dir_results+'tp53_genotype_' + suffix + '_individuals.csv',na_rep='NaN')


# Main

def main():

    two_primer_stats_all = Two_primer_stats_all()

    gt_type = 'wtmut'
    gt_method = 'GL'  # genotype likelihood (GL)

    # Get TP53 genotype data (from scAmp-seq)
    print('Get TP53 genotyoe data (from scAmp-seq)')
    df_tp53_files_wtmut = get_tp53_files(gt_type, gt_method)

    # Get barcodes from adata (from scRNA-seq)
    print('Get barcodes from adata (from scRNA-seq)')
    df_bc = get_bcs_from_adata()

    # Add TP53 genotype data for each variant to adata object
    print('Add TP53 genotype data for each variant to adata object')
    add_gt_to_adata(df_tp53_files_wtmut, gt_type, df_bc, two_primer_stats_all)

    # Add overall TP53 genotype data to adata object
    # Marking cells of E2, SDS1 and SDS2 as WT as they have no TP53 clones
    print('Add overall TP53 genotype data to adata object')
    add_overall_genotype(df_tp53_files_wtmut, gt_type, ['E2','SDS1','SDS2'])

    # Write stats to files
    print('Write stats to files')
    write_overall_stats(gt_type)
    write_overall_stats('wtmut_noTP53clone_wt')
    two_primer_stats_all.write_to_file()

    # Save adata
    adata_file_out = 'adata_filtered_seeds_scANVI_UMAP_unclear-cell-types_cell-cycle_TP53-gt.h5ad'
    adata.write(os.path.join(dir_results, adata_file_out))


if __name__ == "__main__":
    main()
