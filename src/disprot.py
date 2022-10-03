import pandas as pd
import config as cfg
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

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
    # print(set(transformed_regions))
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


def content_count_col_maker(input_df):
    ## because of overlapping regions, we have to count the ranges in a set to discard the repeated numbers
    array_is_in = []
    content_count_array = []
    for index, row in input_df.iterrows():
        set_disorder_region = expand_regions(row.startend)
        content_count_array.append(len(set_disorder_region))  # temp set of data, convert each startend lst to a set,
    input_df['content_count'] = content_count_array
    return input_df


def disprot_prepare_cols():
    dis = pd.read_csv(cfg.data['dis'] + '/DisProt release_2022_06 with_ambiguous_evidences.tsv', sep='\t')
    dis = dis.loc[dis['organism'] == 'Homo sapiens']
    dis['startend'] = dis['start'].astype('string') + '..' + dis['end'].astype('string')
    dis = dis.drop(
        columns=['organism', 'ncbi_taxon_id', 'start', 'end', 'region_sequence', 'confidence', 'obsolete'])
    dis = dis.groupby(['acc', 'name', 'disprot_id', 'term_namespace', 'term'])['startend'].apply(
        ','.join).reset_index()
    dis = content_count_col_maker(dis)
    ## merging with protein lengths from uniprot on 16th sep 2022
    dis_len = pd.read_csv(cfg.data['dis'] + '/protein-lengths-human.tsv', sep='\t')
    dis = pd.merge(dis, dis_len, on='acc')
    dis['content_fraction'] = dis['content_count'] / dis['length']
    return dis


def box_plotter(input_df):
    df = input_df.pivot(index=['index', 'acc'], columns='phenotype', values='content_fraction')

    plt.figure(figsize=(60, 60))  # bigger figsize to have xticklabels shown
    g = sns.catplot(data=df, kind="box").set(title='Disorder content among phenotypes', xlabel='Phenotypes',
                                             ylabel='Disorder content fraction')
    # ylim = (-5, ylim)
    sns.set_style("ticks")
    g.set_xticklabels(rotation=45, va="center", position=(0, -0.02))
    plt.tight_layout()
    plt.savefig(cfg.plots['bp'] + '/Disorder_content_disprot.png', dpi=200)
    plt.show()
    plt.close('all')
    return


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


if __name__ == '__main__':
    hc = pd.read_csv(cfg.data['hc'] + '/smaller-hc-with-phens-column')
    ## Disprot downloaded at sep 14th 2022
    disprot = disprot_prepare_cols()
    dis_hc = pd.merge(disprot, hc, on='acc')
    z = dis_hc.groupby(['term', 'phenotype']).count()  ## to check for the more existing features
    dis_dis = dis_hc.loc[dis_hc['term'] == 'IDPO:00076']
    dis_dis = dis_dis.reset_index()
    dis_dis.to_csv(cfg.data['dis'] + '/disprot_disorder_for_dis_content_boxplot_r.csv') # for Rstudio
    dis_bind = dis_hc.loc[dis_hc['term'] == 'GO:0005515']  # protein binding
    dis_bind = dis_bind.reset_index()
    dis_ord = dis_hc.loc[dis_hc['term'] == 'IDPO:00050']  # disorde to order
    dis_ord = dis_ord.reset_index()
    # box_plotter(dis_dis)
    # sort phenotypes
    # merging with Clinvar to check vars in IDR
    # dis_dis = pd.merge(clinvar_mut_data_maker('patho'), dis_dis, left_on='GeneSymbol', right_on='Gene_name', how='inner')
    dis_bind = pd.merge(clinvar_mut_data_maker('vus'), dis_bind, left_on='GeneSymbol', right_on='Gene_name', how='inner')
    # dis_dis = var_cnt_residue_normaliezer(var_countcol_creator(mutidr_bool_array_maker(dis_dis, 'isin_idr')))
    # dis_dis.to_csv(cfg.data['dis'] + '/Disprot-disorder-variants-in-out-checked.csv')
    dis_bind = var_cnt_residue_normaliezer(var_countcol_creator(mutidr_bool_array_maker(dis_bind, 'isin_idr')))
    dis_bind.to_csv(cfg.data['dis'] + '/Disprot-binding-VUS_variants-in-out-checked.csv')

