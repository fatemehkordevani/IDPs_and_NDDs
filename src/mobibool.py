import sys
import pandas as pd
import numpy as np
from pandas import DataFrame
import config as cfg
import brain as bd
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
import seaborn as sns
import phenotypes as ph


def mobi_phens_col_maker(df1_mobi, df2, df3):
    # mobidb
    # df1_mobi = df1_mobi.drop(columns='start..end')
    subdf1 = df1_mobi[['acc']]
    subdf1['phenotype'] = 'Human'
    # brain
    df2['phenotype'] = 'Brain'
    # ndd
    df3 = df3.drop_duplicates()
    all_phens = subdf1.append([df2, df3]).drop_duplicates()
    new_mobi_all_phens = df1_mobi.merge(all_phens, how='left', on='acc')
    return new_mobi_all_phens


def multidx_df_maker(input_dfs_lst, idx_lst):
    # multi-level index Phenotypes: from: https://www.youtube.com/watch?v=tcRGa2soc-c
    # supposed to do this:
    # mobi_disorder_df = mobi_feature_df.groupby(
    # ['acc', 'feature', 'phenotype']).content_fraction.mean().unstack().sort_index()
    cf_multidx_df = input_dfs_lst[0].groupby(idx_lst).content_fraction.mean().unstack().sort_index()
    cf_multidx_df = cf_multidx_df * 100  # converting cf to percentage
    cc_multidx_df = input_dfs_lst[0].groupby(idx_lst).content_count.mean().unstack().sort_index()
    len_multidx_df = input_dfs_lst[1].groupby(['acc', 'phenotype']).length.mean().unstack().sort_index()  # this is to
    # avoid duplications cuz length is the same and does not differ with mobidb features,
    # so the features should not be taken into account
    return cf_multidx_df, cc_multidx_df, len_multidx_df


def box_plotter(data, title, ylabel, save_route):
    plt.figure(figsize=(60, 60))  # bigger figsize to have xticklabels shown
    g = sns.catplot(data=data, kind="box").set(title=title, xlabel='Phenotypes', ylabel=ylabel)
    # ylim = (-5, ylim)
    sns.set_style("ticks")
    g.set_xticklabels(rotation=45, va="center", position=(0, -0.02))
    plt.tight_layout()
    plt.savefig(save_route)
    plt.close('all')
    return


def violin_plotter(data, title, save_route, ylabel, ylim):
    plt.figure(figsize=(12, 6))
    sns.set_theme(style="whitegrid")
    g = sns.violinplot(data=data, split=True, bw=.1).set(title=title, xlabel='Phenotypes', ylabel=ylabel)
    sns.set_style("ticks")
    # g.set_xticklabels(rotation=45, va="center", position=(0, -0.02))
    plt.ylim(-5, ylim)
    plt.tight_layout()
    plt.savefig(save_route)
    plt.close('all')
    return


def draw_barplot(x, y, data, xticklabel, yscale, save_rout):  # input is DF, not list
    plt.figure(figsize=(12, 6))  # bigger figsize to have xticklabels shown
    sns.set_style("ticks")
    g = sns.barplot(x=x, y=y, data=data)
    # sns.despine(trim=True, offset=2)
    g.set_xticklabels(xticklabel, rotation=0, va="center", position=(0, -0.02))
    sns.color_palette("pastel")
    plt.yscale(yscale)
    plt.tight_layout()
    plt.savefig(save_rout)
    plt.show()
    return


def sig_pep_percent_df_maker():
    # this gets a df with the selected peptides count and then divides them by all proteins of that phenotype or
    # brn, human, ... and creates percentage, all Prs count is taken from length_df cuz it does not have redundancy
    sig_peptide_subdf = pd.read_csv(
        cfg.data['desc-cc'] + '/prediction-signal_peptide-uniprot-cc.csv',
        usecols=['phenotype', 'count'])
    all_phen_pr_count_df = pd.read_csv(cfg.data['desc-len'] + '/length-stats.csv', usecols=['phenotype', 'count'])
    sig_pep_mrged_df = pd.merge(sig_peptide_subdf, all_phen_pr_count_df, on='phenotype')
    sig_pep_mrged_df = sig_pep_mrged_df.rename(
        columns={'phenotype': 'Phenotypes', 'count_x': 'Signal peptide count', 'count_y': 'count_all'})
    sig_pep_mrged_df['Signal peptide percentage'] = (sig_pep_mrged_df['Signal peptide count'] * 100) / \
                                                    sig_pep_mrged_df['count_all']
    return sig_pep_mrged_df


def transmem_pr_percent_df_maker():
    all_phen_pr_count_df = pd.read_csv(cfg.data['desc-len'] + '/length-stats.csv', usecols=['phenotype', 'count'])
    transmemb_subdf = pd.read_csv(
        cfg.data['desc-cc'] + '/prediction-transmembrane-uniprot-cc.csv',
        usecols=['phenotype', 'count'])
    transmem_mrg_df = pd.merge(transmemb_subdf, all_phen_pr_count_df, on='phenotype')
    transmem_mrg_df = transmem_mrg_df.rename(
        columns={'phenotype': 'Phenotypes', 'count_x': 'Transmembrane protein count', 'count_y': 'count_all'})
    transmem_mrg_df['Transmembrane protein percentage'] = (transmem_mrg_df['Transmembrane protein count'] * 100) / \
                                                          transmem_mrg_df['count_all']
    return transmem_mrg_df


def hc_all_phens_generator(df):
    phens_dict = {'ID': ['SysID primary genes (updated 2021-11-18)'], 'Ep': ['Epilepsy ALL genesv2'],
                  'ASD': ['SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)']}
    for k, v in phens_dict.items():
        for i in v:
            df.loc[df[i] == True, i] = k
    # all phens in one column separated by comma
    df['all_phens'] = df[df.columns[3:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    # phens split into multiple rows
    df = (df.set_index(['Unnamed: 0', 'Gene_name', 'acc',
                        'SysID primary genes (updated 2021-11-18)', 'Epilepsy ALL genesv2',
                        'SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)']).apply(
        lambda x: x.str.split(',').explode()).reset_index())
    df = df.loc[df['all_phens'] != 'False'].reset_index()
    df = df.drop(columns=['index', 'Unnamed: 0', 'SysID primary genes (updated 2021-11-18)', 'Epilepsy ALL genesv2',
                          'SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)'])
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


if __name__ == '__main__':
    ## selected features
    features_lst = ['prediction-lip-anchor', 'homology-domain-merge', 'derived-binding_mode_disorder_to_disorder-mobi',
                    'derived-binding_mode_disorder_to_order-mobi', 'prediction-disorder-th_50',
                    'prediction-signal_peptide-uniprot', 'prediction-transmembrane-uniprot']
    titles_lst = ['Linear Interacting Peptides - Anchor', 'Protein domains', 'Binding mode - disorder to disorder',
                  'Binding mode - disorder to order', 'Disorder - majority', 'Signal peptides - Uniprot',
                  'Transmembrane helices']
    titles_lst = [[i] for i in titles_lst]
    # this is a dict with feature names as key and their plot titles as value (then dfs for cc will be added as values)
    feature_dict = dict(zip(features_lst, titles_lst))
    cc_lim_lst = [1010, 810, 1210, 1810, 75, 610, 900, 360, 105, 510, 1810]
    cc_lim_feature_dict = dict(zip(features_lst, cc_lim_lst))
    cf_lim_lst = [105, 65, 105, 105, 105, 105, 105, 105, 105, 105, 105]
    # cf_lim_lst = [110, 65, 110, 110, 110, 110, 110, 110, 110, 110, 110] for violins I used 110 instead of 105
    cf_lim_feature_dict = dict(zip(features_lst, cf_lim_lst))
    # now add content_count limit of each feature to your dict as second value (idx=1)
    # and content fraction limit as idx=2
    for feature in features_lst:
        feature_dict[feature].append(cc_lim_feature_dict[feature])
        feature_dict[feature].append(cf_lim_feature_dict[feature])
    # in the end use decorator thing with the @
    phens_lst = ['Human', 'Brain', 'ASD', 'Ep', 'ID', 'ADHD', 'SCZ', 'Cancer', 'T2D']
    ## import dfs # (mobidb)
    mobidb = pd.read_csv(cfg.data['clin'] + '/mobidb_result_2022-08-08T13_12_39.379Z.tsv', sep='\t')
    ## NDD
    phens_subdf = pd.read_csv(cfg.data['hc'] + '/smaller-hc-with-phens-column')
    ndd_pr_lst = phens_subdf['acc'].unique().tolist()  # 1597
    ## brain
    brain_prot_lst = bd.brain_pr_lst_generator()  # n: 8320
    brain_subdf = DataFrame(brain_prot_lst, columns=['acc'])
    # mutual_brain_ndd_prs_lst = [i for i in brain_prot_lst if i in ndd_pr_lst]  # 455
    # # ## Mobidb, Brain and Ndd dfs with all variations (in/out of IDR)
    # # mobidb, brain_subdf, phens_subdf = var.all_vars_or_vars_inidr('all')
    # #
    # new mobidb with one column for phens of NDD, human, and brain
    mobidb = mobi_phens_col_maker(mobidb, brain_subdf, phens_subdf)
    mobidb.to_csv(cfg.data['vars'] + '/all-vars-mobidb-plus-phenotype-column.csv')
    #
    mobi_feature_df = mobidb[mobidb.feature.isin(features_lst)]
    # # multi-idx-dfs
    mobi_disorder_df, mobi_cont_count_df, mobi_length_df = multidx_df_maker(
        [mobi_feature_df, mobidb], ['acc', 'feature', 'phenotype'])

    ## dfs for statistical tests
    disorder_stat = mobi_disorder_df.loc[(slice(None), 'prediction-disorder-th_50'), phens_lst]
    disorder_stat.to_csv(cfg.data['stats'] + '/mobidb-disorder-for-stats')
    lip_stat = disorder_stat = mobi_disorder_df.loc[(slice(None), 'prediction-lip-anchor'), phens_lst]
    lip_stat.to_csv(cfg.data['stats'] + '/mobidb-lip-for-stats')
    domain_stat = disorder_stat = mobi_disorder_df.loc[(slice(None), 'homology-domain-merge'), phens_lst]
    domain_stat.to_csv(cfg.data['stats'] + '/mobidb-domain-for-stats')

    ## Boxplots
    # content count
    for key in feature_dict.keys():
        box_plotter(data=mobi_cont_count_df.loc[(slice(None), key), phens_lst],
                    save_route=(cfg.plots['bp-cc'] + '/' + key + '-cc' + '1.png'),
                    title=feature_dict[key][0], ylabel='Residues count')
    # content fraction
    for key in feature_dict.keys():
        box_plotter(data=mobi_disorder_df.loc[(slice(None), key), phens_lst],
                    save_route=(cfg.plots['bp-cf'] + '/' + key + '-cf' + '1.png'),
                    title=feature_dict[key][0], ylabel='Content (%)')
    # Length
    box_plotter(data=mobi_length_df.loc[(slice(None)), phens_lst],
                save_route=(cfg.plots['bp-len'] + '/' + 'len4200' + '1.png'),
                title='Protein sequence length', ylabel='Residues count')
    ## Violin plots
    # content count
    for key in feature_dict.keys():
        violin_plotter(data=mobi_cont_count_df.loc[(slice(None), key), phens_lst],
                       save_route=(cfg.plots['vp-cc'] + '/' + key + '-cc' + '1.png'),
                       title=feature_dict[key][0], ylabel='Residues count', ylim=feature_dict[key][1])
    # content fraction
    for key in feature_dict.keys():
        violin_plotter(data=mobi_disorder_df.loc[(slice(None), key), phens_lst],
                       save_route=(cfg.plots['vp-cf'] + '/' + key + '-cf' + '1.png'),
                       title=feature_dict[key][0], ylabel='Content (%)', ylim=feature_dict[key][2])
    # Length
    violin_plotter(data=mobi_length_df.loc[(slice(None)), phens_lst],
                   save_route=(cfg.plots['vp-len'] + '/' + 'len4200' + '1.png'),
                   title='Protein sequence length', ylabel='Residues count', ylim=4200)

    ## writing data statistics to CSV
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    for each_f in features_lst:
        mobi_disorder_df.loc[(slice(None), each_f), phens_lst].describe().T. \
            to_csv(cfg.data['desc-cf'] + '/' + each_f + '-cf.csv')
        mobi_cont_count_df.loc[(slice(None), each_f), phens_lst].describe().T. \
            to_csv(cfg.data['desc-cc'] + '/' + each_f + '-cc.csv')
    mobi_length_df.loc[slice(None), phens_lst].describe().T.to_csv(cfg.data['desc-len'] + '/length-stats.csv')

    sig_pep_percent_df = sig_pep_percent_df_maker()
    transmem_pr_percent_df = transmem_pr_percent_df_maker()

    ## Barplots
    # Signal peptide
    draw_barplot(x='Phenotypes', y='Signal peptide percentage', data=sig_pep_percent_df, xticklabel=phens_lst,
                 yscale='linear', save_rout=cfg.plots['bar-sptm'] + '/sig-peptide-percent.png')
    draw_barplot(x='Phenotypes', y='Signal peptide count', data=sig_pep_percent_df, xticklabel=phens_lst,
                 yscale='log', save_rout=cfg.plots['bar-sptm'] + '/sig-peptide-count.png')
    # Transmembrane protein
    draw_barplot(x='Phenotypes', y='Transmembrane protein percentage', data=transmem_pr_percent_df,
                 xticklabel=phens_lst,
                 yscale='linear', save_rout=cfg.plots['bar-sptm'] + '/transmemb-prots-percent.png')
    draw_barplot(x='Phenotypes', y='Transmembrane protein count', data=transmem_pr_percent_df, xticklabel=phens_lst,
                 yscale='log', save_rout=cfg.plots['bar-sptm'] + '/transmemb-prots-count.png')
