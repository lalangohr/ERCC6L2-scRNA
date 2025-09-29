"""
Part 8 of scRNA-seq data analysis. 

Pathway analysis of DE genes identified in part 7 of scRNA-seq data analysis for cell type annotation of unclear cell types.

Adapted from
https://maayanlab.cloud/Enrichr/help#api
https://github.com/MaayanLab/enrichr_api/blob/master/main.py
https://github.com/wdecoster/enrichr_cli/blob/master/enrichr-cli.py
"""

import json
import requests
import os
import csv
import pandas as pd
from time import sleep
from pathlib import Path


# Folder in which result is written
dir_results = './results'


def main():

    files_in = ['wilcoxon_celltype_unclear_MkP_cluster0_using_subset.csv',
                'wilcoxon_celltype_unclear_MkP_cluster1_using_subset.csv',
                'wilcoxon_celltype_unclear_myelo-mono_cluster3_using_subset.csv']

    library = 'Azimuth_2023'

    for file_in in files_in:

        file_out = 'enrichR_' + library + '_' + file_in

        df_de = pd.read_csv(os.path.join(dir_results,file_in))

        # Use only upregulated genes with log2FC>=1.0 and p_adj<=0.05:
        degs = df_de[(df_de['pvals_adj']<=0.05) & (df_de['logfoldchanges']>=1.0)].names

        # Analyze gene set. Returns JSON object with unique ID for analysis results.

        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'

        genes_str = '\n'.join(degs)
        description = 'Gene list'
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }

        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        data_id = json.loads(response.text)
        print(data_id)

        # Download file of enrichment results. Returns text file of enrichment analysis results.

        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/export'

        query_string = '?userListId=%s&filename=%s&backgroundType=%s'
        filename = 'example_enrichment'
        user_list_id = data_id['userListId']

        gene_set_library = library

        url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)
        response = requests.get(url, stream=True)

        with open(os.path.join(dir_results,file_out), 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024): 
                if chunk:
                    f.write(chunk)


if __name__ == "__main__":
    main()
