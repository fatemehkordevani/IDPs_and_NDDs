import pandas as pd
import config as cfg
import numpy as np


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


def clinvar_mut_data_maker(clin_sig_type):
    ## clinvar version: variant summary.txt 135,193,289	2022-08-01 17:10:07	2022-08-01 17:10:07

    clinvar = pd.read_csv(cfg.data['clin'] + '/variant_summary.txt', sep='\t', low_memory=False)
    clinvar = clinvar.drop(
        columns=['Assembly', 'ChromosomeAccession', 'Start', 'Stop', 'PositionVCF', 'ReferenceAlleleVCF',
                 'AlternateAlleleVCF'])
    clinvar = clinvar.drop_duplicates()
    if clin_sig_type == 'patho':
        clin_sig = ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic', 'Pathogenic; association',
                    'Likely pathogenic; risk factor', 'Pathogenic; risk factor', 'Pathogenic; Affects',
                    'Pathogenic/Likely pathogenic; risk factor', 'Likely pathogenic; Affects',
                    'Pathogenic; drug response', 'Likely pathogenic; drug response', 'Pathogenic; protective',
                    'Pathogenic/Likely pathogenic; drug response', 'Pathogenic; association; protective',
                    'Pathogenic; confers sensitivity', 'Likely pathogenic; association']
    elif clin_sig_type == 'vus':
        clin_sig = ['Uncertain significance', 'Uncertain significance; risk factor',
                    'Uncertain significance; drug response', 'Uncertain significance; Pathogenic/Likely pathogenic',
                    'Uncertain significance; other', 'Uncertain significance; association',
                    'Uncertain significance; Affects']
    clinvar = clinvar.loc[clinvar.ClinicalSignificance.isin(clin_sig)]
    # clinvar = clinvar.loc[(clinvar['ClinicalSignificance'] == 'Pathogenic') |
    #                       (clinvar['ClinicalSignificance'] == 'Likely pathogenic')
    #                       | (clinvar['ClinicalSignificance'] == 'Pathogenic/Likely pathogenic')]
    # clinvar = clinvar.loc[clinvar['PhenotypeList'] != 'not provided']
    ## getting mutation positions
    # clinvar[['Name', 'pr_change']] = clinvar['Name'].str.split('\(p\.', 1, expand=True)
    # clinvar = clinvar.dropna(subset=['pr_change'])
    clinvar['pr_change'] = clinvar.Name.str.extract('(p\.([A-Z])([a-z])([a-z])(\d+)(.+))', expand=False)[0]
    # clinvar['pr_change'] = clinvar['pr_change'].str.split('p.', 1, expand=True)[0]
    clinvar['pr_change'] = clinvar.pr_change.str.replace(')', '')
    clinvar['pr_change'] = clinvar.pr_change.str.replace('p\.', '')
    clinvar = clinvar.dropna(subset=['pr_change'])
    # clinvar = clinvar.loc[~clinvar.pr_change.str.contains('=')]
    clinvar['ref_aa'] = clinvar['pr_change'].str[:3]
    clinvar['alt_aa'] = clinvar['pr_change'].str[-3:]
    clinvar['pr_pos'] = clinvar.pr_change.str.extract('(\d+)')
    clinvar = clinvar.dropna(subset=['pr_pos'])
    clinvar = clinvar.drop_duplicates()
    clinvar['var_id'] = clinvar.index + 1000
    clinvar['var_id'] = 'var_' + clinvar['var_id'].astype(str)
    return clinvar


def var_in_idr_lip_dfs(clin_phen, mobi):
    mobi_dis = mobi.loc[mobi['feature'] == 'prediction-disorder-th_50']  # (77629, 6)
    clin_mobi_ndd = pd.merge(clin_phen, mobi_dis, on='acc', how='inner')  # (35260, 37)
    ## check if modification position is in idr
    clin_mobi_ndd = mutidr_bool_array_maker(clin_mobi_ndd, 'isin_idr')
    ## count number of vars and var/dis ratio
    clin_mobi_ndd = var_cnt_residue_normaliezer(var_countcol_creator(clin_mobi_ndd))
    clin_mobi_ndd.to_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb-VUS.csv')
    in_idr_vars = clin_mobi_ndd.loc[clin_mobi_ndd['isin_idr'] == 1]
    # clin_mobidb_ndd = pd.read_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb.csv')
    in_idr_var_ids = in_idr_vars['var_id'].unique().tolist()
    ## Vars in LIPs
    mobi_lip = mobi.loc[mobi['feature'] == 'prediction-lip-anchor']
    ndd_mobi_lip = pd.merge(clin_phen, mobi_lip, left_on='acc', right_on='acc', how='inner')  # ()
    ndd_mobi_lip = ndd_mobi_lip.loc[ndd_mobi_lip.var_id.isin(in_idr_var_ids)]
    ndd_mobi_lip = mutidr_bool_array_maker(ndd_mobi_lip, 'isin_lip')
    ndd_mobi_lip.to_csv(cfg.data['hc'] + '/mobidb_vars_in_lips_checked-VUS.csv')
    return clin_mobi_ndd, ndd_mobi_lip


if __name__ == '__main__':
    hc = pd.read_csv(cfg.data['hc'] + '/smaller-hc-with-phens-column')
    clinvar = clinvar_mut_data_maker('patho')
    clinvus = clinvar_mut_data_maker('vus')
    ## merge my protein list with clinvar
    # clinvar
    clinvar_ndd = pd.merge(hc, clinvar, left_on='Gene_name', right_on='GeneSymbol')
    clinvar_ndd.to_csv(cfg.data['hc'] + '/clinvar_and_phenotypes_mutation_positions.csv')
    clinvar_ndd = pd.read_csv(cfg.data['hc'] + '/clinvar_and_phenotypes_mutation_positions.csv')
    clinvar_ndd = clinvar_ndd.drop(columns=['Unnamed: 0.1', 'Unnamed: 0'])
    # clinvus
    clinvus_ndd = pd.merge(hc, clinvus, left_on='Gene_name', right_on='GeneSymbol')
    clinvus_ndd.to_csv(cfg.data['hc'] + '/clinvus_and_phenotypes_mutation_positions.csv')
    clinvus_ndd = pd.read_csv(cfg.data['hc'] + '/clinvus_and_phenotypes_mutation_positions.csv')
    clinvus_ndd = clinvus_ndd.drop(columns=['Unnamed: 0.1', 'Unnamed: 0'])
    ## Mobidb downloaded on Aug 8th 2022, 78,106 entries (but 77629 prs after filtering for disorder consensus)
    mobidb = pd.read_csv(cfg.data['clin'] + '/mobidb_result_2022-08-08T13_12_39.379Z.tsv', sep='\t')
    mobidb = mobidb.rename(columns={'start..end': 'startend'})

    ## merge clinvar data positions and stuff with mobidb disorder
    # clin_mobidb_ndd, ndd_mobidb_lip = var_in_idr_lip_dfs(clinvar_ndd, mobidb)
    clin_mobidb_ndd = pd.read_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb.csv')
    ndd_mobidb_lip = pd.read_csv(cfg.data['hc'] + '/mobidb_vars_in_lips_checked.csv')
    in_idr_vars = clin_mobidb_ndd.loc[clin_mobidb_ndd['isin_idr'] == 1]
    # clinvus
    vus_mobidb_ndd, vus_mobidb_lip = var_in_idr_lip_dfs(clinvus_ndd, mobidb)
    vus_mobidb_ndd = pd.read_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb-VUS.csv')
    vus_mobidb_lip = pd.read_csv(cfg.data['hc'] + '/mobidb_vars_in_lips_checked-VUS.csv')
    ## var in ptm
    # var_ptm_checked_df = dismaj_var_in_ptm_df_generator(in_idr_vars)
    # var_ptm_checked_df.to_csv(cfg.data['ptm'] + '/in_idr_vars_ptm_checked.csv')
    # var_ptm_checked_df = pd.read_csv(cfg.data['ptm'] + '/in_idr_vars_ptm_checked.csv')
    # ex = var_ptm_checked_df.loc[(var_ptm_checked_df['ptm_type'] == 'disulfide bond') & (var_ptm_checked_df['var_in_ptm']== 1)]
    #
    # cancer_pr = clinvar_ndd.loc[clinvar_ndd['phenotype'] == 'Cancer', 'acc'].unique().tolist()
    # other_pr = clinvar_ndd.loc[clinvar_ndd['phenotype'] != 'Cancer', 'acc'].unique().tolist()
