# check brain proteome data


import pandas as pd
import config as cfg
import sys
import re


def delimiter_split(delimiter, list):
    new_lst = []
    for i in list:
        if delimiter in i:
            new_lst.extend(i.split(delimiter))
        elif delimiter not in i:
            new_lst.append(i)
    return new_lst


def brain_pr_lst_generator():
    # Data from the Human Protein Atlas in tab-separated format
    # This file contains a subset of the data in the Human Protein Atlas version 20.1
    brain_df = pd.read_table(cfg.data['brain'] + '/proteinatlas.tsv', error_bad_lines=False, sep='\t')

    # Normal tissue data
    #  The data is based on The Human Protein Atlas version 20.1 and Ensembl version 92.38.
    # (1118517, 6)
    normal_tissue_df = pd.read_table(cfg.data['brain'] + '/normal_tissue.tsv', error_bad_lines=False, sep='\t')
    norm_tissue_cols_lst = normal_tissue_df.columns.tolist()

    brain_areas_lst = ['hypothalamus', 'hippocampus', 'cerebral cortex', 'cerebellum', 'choroid plexus',
                       'dorsal raphe', 'pituitary gland', 'caudate']
    # to be filtered based on only brain parts (discarding other tissues)
    brain_filtered_df = pd.DataFrame(columns=norm_tissue_cols_lst)  # defining empty df to be filled
    for each_area in brain_areas_lst:
        df_temp = normal_tissue_df[(normal_tissue_df['Tissue'] == each_area)]
        brain_filtered_df = brain_filtered_df.append(df_temp)  # (161558, 6)

    # filter by level of expression  (74451, 6)
    brain_filtered_df1 = brain_filtered_df[(brain_filtered_df['Level'] == 'Low') |
                                           (brain_filtered_df['Level'] == 'Medium') |
                                           (brain_filtered_df['Level'] == 'High')]
    # filter by Reliability  (60461, 6)
    brain_filtered_df1 = brain_filtered_df1[(brain_filtered_df['Reliability'] == 'Approved') |
                                            (brain_filtered_df['Reliability'] == 'Enhanced') |
                                            (brain_filtered_df['Reliability'] == 'Supported')]
    ## this df contains Gene names and Ensembl IDs which I'll get the unique ones as list and IDmap
    # from uniprot for IDs
    gene_ensembl_lst = brain_filtered_df1['Gene'].unique()
    # since the number of ensembl IDs are not that much, let's check them first in the first df here (brain df)
    subdf_toget_acc = brain_df[brain_df.Ensembl.isin(gene_ensembl_lst)]
    brain_uniprot_lst = subdf_toget_acc['Uniprot'].tolist()
    # drop redundancies
    brain_uniprot_lst = list(set(brain_uniprot_lst))
    brain_uniprot_lst = [x for x in brain_uniprot_lst if str(x) != 'nan']
    brain_uniprot_lst = delimiter_split(', ', brain_uniprot_lst)

    return brain_uniprot_lst


if __name__ == '__main__':
    # if this number is not enough, you can try to include other brain dataframes as well, e.g: Brainspan
    # brain_pr_lst = brain_pr_lst_generator()
    brain_df = pd.read_table(cfg.data['brain'] + '/proteinatlas.tsv', error_bad_lines=False, sep='\t')
    normal_tissue_df = pd.read_table(cfg.data['brain'] + '/normal_tissue.tsv', error_bad_lines=False, sep='\t')
    normal_tissue_lst = normal_tissue_df['Level'].unique().tolist()
    br_new = pd.read_table(cfg.data['brain'] + '/brain_category_rna_Any_Region.tsv', sep='\t')
