import pandas as pd
import config as cfg
import numpy as np
from fractions import Fraction

# downloaded from https://genetrek.pasteur.fr/ on August 23rd 2022
# hc_genes = pd.read_csv(cfg.data['hc'] + '/genetrek-data-v6-2022-03-31.tsv', sep='\t')
# hc_genes = hc_genes.loc[hc_genes['curatedLists.highConfidenceNddV2'] == True]

def expand_regions(region_ranges_lst):
    transformed_regions = []
    region_ranges_lst = region_ranges_lst.split(',')
    # region_ranges_lst = list(region_ranges_lst)
    for region in region_ranges_lst:
        start = int(region.split('..')[0])
        end = int(region.split('..')[1])
        while start <= end:
            transformed_regions.append(start)
            start += 1
    return set(transformed_regions)


def mutidr_bool_array_maker(input_df, checked_col_name):
    # input df is merged df of mobidb and mutation positions from mutinfo
    ## checks if mutation position is in startend disorder region of mobidb or not
    array_is_in = []  # will be filled with boolean of 0,1
    for index, row in input_df.iterrows():
        set_disorder_region = expand_regions(row.startend)  # temp set of data, convert each startend lst to a set,

        if int(row.pr_pos) in set_disorder_region:
            array_is_in.append('1')
        else:
            array_is_in.append('0')
    input_df[checked_col_name] = array_is_in
    return input_df


## related to vars in PTM
def uniprot_var_in_ptm_checker(df):  # input is merged 1-variants+in/out disorder columns and 2-ptms df from uniprot
    array_isin = []
    for index, row in df.iterrows():
        if '&' not in row.ptm_pos:
            # here we fist make a list of ptm_pos with range of 3 AAs before and after the ptm-point (not for disulfide)
            ptm_pos_tmp_lst = []
            ptm_pos_tmp_lst = np.arange(int(row.ptm_pos) - 3, int(row.ptm_pos) + 4)
            # (delete negative numbers later)
            if int(row.pr_pos) in ptm_pos_tmp_lst:
                array_isin.append('1')
            else:
                array_isin.append('0')
        elif '&' in row.ptm_pos:
            set_ptm_pos = set(str(row.ptm_pos).split('&'))
            if str(row.pr_pos) in set_ptm_pos:
                array_isin.append('1')
            else:
                array_isin.append('0')
    df['var_in_ptm'] = array_isin
    return df


def var_in_ptm_checker(input_df, source):
    if source == 'uniprot':
        return uniprot_var_in_ptm_checker(input_df)


def dismaj_var_in_ptm_df_generator(input_df_mobi_var_in):
    # this is just to organize the code, disorder_majority and ptms df are being merged and then with
    # var_in_ptm_checker a new df is produce that shows var positions that are in ptm sites (ptm pos)
    ptms_df = pd.read_csv(cfg.data['ptm-u'] + '/uniprot-ptms-all.csv',
                          usecols=['acc', 'ptm_pos', 'description', 'ptm_type'])  # (140538, 4)
    dismaj_ptm_df = pd.merge(input_df_mobi_var_in, ptms_df, on='acc')
    ptm_checked_dismaj_df = var_in_ptm_checker(dismaj_ptm_df, 'uniprot')
    ptm_checked_dismaj_df.to_csv(cfg.data['ptm'] + '/uniprot-vars-inptm-checked.csv')
    return ptm_checked_dismaj_df


def var_countcol_creator(df):  # a percentage column for variants in IDR / all Vars
    # input it the out put of df_feature_filtered method
    df['isin_idr'] = df['isin_idr'].apply(pd.to_numeric)
    var_count_dic = (df.pivot_table(columns=['acc'], aggfunc='size')).to_dict()
    # dict values to list so that I could append other data to dict values
    for key, value in var_count_dic.items():
        var_count_dic[key] = [value]
    for k in var_count_dic:  # adds the count of each acc var being in idr
        each_acc_vars_in_idr_count = df.loc[df['acc'] == k, 'isin_idr'].sum()
        var_count_dic[k].append(each_acc_vars_in_idr_count)
    # dict to df
    var_count_df = pd.DataFrame.from_dict(var_count_dic, orient='index', columns=['total_vars', 'in_idr_vars'])
    var_count_df = var_count_df.reset_index()
    var_count_df['out_idr_vars'] = (var_count_df['total_vars'] - var_count_df['in_idr_vars'])
    var_count_df = var_count_df.reset_index()
    var_count_df = var_count_df.rename(columns={'index': 'acc'})
    mrg_var_and_count_df = pd.merge(df, var_count_df, on='acc')
    # mrg_var_and_count_df = mrg_var_and_count_df.drop(columns=['Unnamed: 0', 'Unnamed: 0_x', 'Unnamed: 0_y'])
    return mrg_var_and_count_df


def var_cnt_residue_normaliezer(df):  # gets df with count columns for vars in/out idr and normalize them based on
    # number of disordered residues or non disordered, respectively
    ## this where method conditions is to prevent nan values for fully disordered proteins (length-content_count == 0)
    # in other proteins cases it just divides number of vars on number of residues
    df['in_idr_vars_perc'] = np.where((df['length'] == df['content_count']), 1,
                                      (df['in_idr_vars'] / df['content_count']))
    df['out_idr_vars_perc'] = np.where((df['length'] == df['content_count']), 0,
                                       (df['out_idr_vars'] / (df['length'] - df['content_count'])))

    # df['in_idr_vars_perc'] = (df['in_idr_vars'] / df['content_count'])
    # df['out_idr_vars_perc'] = (df['out_idr_vars'] / (df['length'] - df['content_count']))
    ## now we will divide each by sum of the calculated data for in+out IDRs to produce complementary % values for cols
    sum_of_normalized_vars = df['in_idr_vars_perc'] + df['out_idr_vars_perc']
    df['in_idr_vars_perc'] = df['in_idr_vars_perc'] / sum_of_normalized_vars
    df['out_idr_vars_perc'] = df['out_idr_vars_perc'] / sum_of_normalized_vars
    return df


def hc_all_phens_generator(df):
    phens_dict = {'ID': ['SysID primary genes (updated 2021-11-18)'], 'Ep': ['Epilepsy ALL genesv2'],
                  'ASD': ['SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)']}
    for k, v in phens_dict.items():
        for i in v:
            df.loc[hc[i] == True, i] = k
    # all phens in one column separated by comma
    df['all_phens'] = df[df.columns[3:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    # phens split into multiple rows
    df = (df.set_index(['Unnamed: 0', 'Gene_name', 'acc',
                        'SysID primary genes (updated 2021-11-18)', 'Epilepsy ALL genesv2',
                        'SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)']).apply(
        lambda x: x.str.split(',').explode()).reset_index())
    df = df.loc[df['all_phens'] != 'False'].reset_index()
    df = df.drop(columns=['index', 'Unnamed: 0', 'SysID primary genes (updated 2021-11-18)', 'Epilepsy ALL genesv2',
                          'SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)', 'Gene_name'])
    df = df.rename(columns={'Protein_acc': 'acc', 'all_phens': 'phenotype'})
    return df


def cosmic_phen_df_maker():
    ## this gives out genes that correspond to only one protein, gets 733 genes and give out 700
    cosmic = pd.read_csv(cfg.data['cos'] + '/Census_allSun Aug 28 14_40_13 2022.tsv', sep='\t')
    cos_genes = cosmic['Gene Symbol'].unique().tolist()
    cos_proteins = pd.read_csv(cfg.data[
                                   'cos'] + '/uniprot-download_true_fields_accession_2Creviewed_2Cprotein_existenc-2022.08.28-14.02.40.32.tsv',
                               sep='\t')
    cos_proteins = cos_proteins.drop(columns=['Protein existence', 'Protein names'])
    cos_proteins = cos_proteins.loc[cos_proteins['Reviewed'] == 'reviewed']
    cos_proteins_count = cos_proteins.groupby('From').count()
    cos_proteins_count = cos_proteins_count.loc[cos_proteins_count['Entry'] == 1]
    unique_gene_protein = cos_proteins_count.reset_index()
    unique_gene_protein = unique_gene_protein['From'].unique().tolist()
    cos_proteins = cos_proteins.loc[cos_proteins.From.isin(unique_gene_protein)]
    del cos_proteins['Reviewed']
    cos_proteins['phenotype'] = 'Cancer'
    cos_proteins = cos_proteins.rename(columns={'From': 'Gene_name', 'Entry': 'acc'})
    return cos_proteins


def clinvar_mut_data_maker():
    ## clinvar version: variant summary.txt 135,193,289	2022-08-01 17:10:07	2022-08-01 17:10:07

    clinvar = pd.read_csv(cfg.data['clin'] + '/variant_summary.txt', sep='\t', low_memory=False)
    clinvar = clinvar.drop(columns=['Assembly', 'ChromosomeAccession', 'Start', 'Stop', 'PositionVCF', 'ReferenceAlleleVCF',
                                    'AlternateAlleleVCF'])
    clinvar = clinvar.drop_duplicates()
    clinvar = clinvar.loc[(clinvar['ClinicalSignificance'] == 'Pathogenic') |
                          (clinvar['ClinicalSignificance'] == 'Likely pathogenic')
                          | (clinvar['ClinicalSignificance'] == 'Pathogenic/Likely pathogenic')]
    clinvar = clinvar.loc[clinvar['PhenotypeList'] != 'not provided']
    clinvar['var_id'] = clinvar.index + 1000
    clinvar['var_id'] = 'var_' + clinvar['var_id'].astype(str)
    print(clinvar)
    ## getting mutation positions
    clinvar[['Name', 'pr_change']] = clinvar['Name'].str.split('p\.', 1, expand=True)
    clinvar = clinvar.dropna(subset=['pr_change'])
    clinvar['pr_change'] = clinvar.pr_change.str.extract('(^([A-Z])([a-z])([a-z])(\d+)(.+))', expand=False)[0]
    # clinvar['pr_change'] = clinvar['pr_change'].str.split('p.', 1, expand=True)[0]
    clinvar['pr_change'] = clinvar.pr_change.str.replace(')', '')
    clinvar = clinvar.dropna(subset=['pr_change'])
    clinvar = clinvar.loc[~clinvar.pr_change.str.contains('=')]
    clinvar['ref_aa'] = clinvar['pr_change'].str[:3]
    clinvar['alt_aa'] = clinvar['pr_change'].str[-3:]
    clinvar['pr_pos'] = clinvar.pr_change.str.extract('(\d+)')
    clinvar = clinvar.dropna(subset=['pr_pos'])
    return clinvar


if __name__ == '__main__':

    hc = pd.read_csv(cfg.data['hc'] + '/genetrek-2022-08-23.tsv', sep='\t')
    gene_protein = pd.read_excel(cfg.data['hc'] + '/HC_NDD_genes list_1737to1721.xlsx')
    ## merge with protein names
    hc = pd.merge(gene_protein, hc, left_on='Gene_name', right_on='Gene', how='inner')
    hc = hc[['Gene_name', 'acc', 'SysID primary genes (updated 2021-11-18)', 'Epilepsy ALL genesv2',
             'SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)']].drop_duplicates()
    hc.to_csv(cfg.data['hc'] + '/smaller-hc-columns.csv')
    hc = pd.read_csv(cfg.data['hc'] + '/smaller-hc-columns.csv')
    ## cosmic genes downloaded August 28th
    cosmic = cosmic_phen_df_maker()
    ##
    hc_phens = hc_all_phens_generator(hc)
    hc = pd.merge(hc, hc_phens, on='acc')
    hc = hc.drop(columns=['Unnamed: 0', 'SysID primary genes (updated 2021-11-18)', 'Epilepsy ALL genesv2',
       'SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)', 'all_phens'])
    ## concatenating hc and cosminc
    hc = pd.concat([hc, cosmic])
    hc.to_csv(cfg.data['hc']+ 'smaller-hc-with-phens-column')
    clinvar = clinvar_mut_data_maker()
    ## merge my protein list with clinvar
    clinvar_ndd = pd.merge(hc, clinvar, left_on='Gene_name', right_on='GeneSymbol')
    clinvar_ndd.to_csv(cfg.data['hc'] + '/clinvar_and_phenotypes_mutation_positions.csv')
    # clinvar_ndd = pd.read_csv(cfg.data['hc'] + '/clinvar_and_phenotypes_mutation_positions.csv')


    ## Mobidb downloaded on Aug 8th 2022, 78,106 entries (but 77629 prs after filtering for disorder consensus)
    mobidb = pd.read_csv(cfg.data['clin'] + '/mobidb_result_2022-08-08T13_12_39.379Z.tsv', sep='\t')
    mobidb = mobidb.rename(columns={'start..end': 'startend'})
    mobidb_disorder = mobidb.loc[mobidb['feature'] == 'prediction-disorder-th_50']  # (77629, 6)
    ## merge clinvar data positions and stuff with mobidb disorder
    clin_mobidb_ndd = pd.merge(clinvar_ndd, mobidb_disorder, on='acc', how='inner')  # (35260, 37)
    ## check if modification position is in idr
    clin_mobidb_ndd = mutidr_bool_array_maker(clin_mobidb_ndd, 'isin_idr')
    ## count number of vars and var/dis ratio
    clin_mobidb_ndd = var_cnt_residue_normaliezer(var_countcol_creator(clin_mobidb_ndd))
    ##
    # clin_mobidb_ndd.to_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb.csv')
    clin_mobidb_ndd = pd.read_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb.csv')
    in_idr_vars = clin_mobidb_ndd.loc[clin_mobidb_ndd['isin_idr'] == 1]
    ## Vars in LIPs
    mobidb_lip = mobidb.loc[mobidb['feature'] == 'prediction-lip-anchor']
    ndd_mobidb_lip = pd.merge(clinvar_ndd, mobidb_lip, left_on='acc', right_on='acc', how='inner')  # ()
    ndd_mobidb_lip = mutidr_bool_array_maker(ndd_mobidb_lip, 'isin_lip')
    ndd_mobidb_lip.to_csv(cfg.data['hc'] + '/mobidb_vars_in_lips_checked.csv')
    in_lip_vars = ndd_mobidb_lip.loc[ndd_mobidb_lip['isin_lip'] == '1']
    ## var in ptm
    var_ptm_checked_df = dismaj_var_in_ptm_df_generator(in_idr_vars)
    var_ptm_checked_df.to_csv(cfg.data['ptm'] + '/in_idr_vars_ptm_checked.csv')
    var_ptm_checked_df = pd.read_csv(cfg.data['ptm'] + '/in_idr_vars_ptm_checked.csv')
    ex = var_ptm_checked_df.loc[(var_ptm_checked_df['ptm_type'] == 'disulfide bond') & (var_ptm_checked_df['var_in_ptm']== 1)]

    cancer_pr = clinvar_ndd.loc[clinvar_ndd['phenotype'] == 'Cancer', 'acc'].unique().tolist()
    other_pr = clinvar_ndd.loc[clinvar_ndd['phenotype'] != 'Cancer', 'acc'].unique().tolist()