import pandas as pd
import config as cfg


def hc_all_phens_generator(df):
    # adding adhd and SCZ to the existing df

    ##
    phens_dict = {'Ep': ['Epilepsy ALL genesv2'], 'ID': ['SysID primary genes (updated 2021-11-18)'],
                  'ASD': ['SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)']}
    for k, v in phens_dict.items():
        for i in v:
            df.loc[df[i] == True, i] = k
    # all phens in one column separated by comma
    df['all_phens'] = df[df.columns[2:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    # phens split into multiple rows
    df = (df.set_index(['Gene_name', 'acc',
                        'SysID primary genes (updated 2021-11-18)', 'Epilepsy ALL genesv2',
                        'SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)']).apply(
        lambda x: x.str.split(',').explode()).reset_index())
    df = df.loc[df['all_phens'] != 'False'].reset_index()
    df = df.drop(columns=['index', 'SysID primary genes (updated 2021-11-18)', 'Epilepsy ALL genesv2',
                          'SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)', 'Gene_name'])
    df = df.rename(columns={'Protein_acc': 'acc', 'all_phens': 'phenotype'})
    return df


def phen_df_maker():
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
    ## import adhd and scz
    adhd_scz = pd.read_csv(cfg.data['ph'] + '/DBD-Genes-Full-Data.csv')
    ## ADHD
    adhd = adhd_scz.loc[adhd_scz['ADHD'] == 'X']
    adhd = adhd['Gene'].unique().tolist()
    adhd = pd.read_csv(cfg.data['ph'] + '/adhd-acc.tsv', sep='\t')
    adhd['phenotype'] = 'ADHD'
    adhd = adhd.rename(columns={'From': 'Gene_name', 'Entry': 'acc'})
    ## SCZ
    scz = adhd_scz.loc[adhd_scz['Schizophrenia'] == 'X']
    scz = scz['Gene'].unique().tolist()
    scz = pd.read_csv(cfg.data['ph'] + '/scz-acc.tsv', sep='\t')
    scz['phenotype'] = 'SCZ'
    scz = scz.rename(columns={'From': 'Gene_name', 'Entry': 'acc'})
    ## Type 2 diabetes, scores 5 and 4
    t2d = pd.read_csv(cfg.data['ph'] + '/t2d-7thSep2022.txt')
    t2d = t2d['Gene'].unique().tolist()
    t2d = pd.read_csv(cfg.data['ph'] + '/t2d-acc.tsv', sep='\t')
    t2d['phenotype'] = 'T2D'
    t2d = t2d.rename(columns={'From': 'Gene_name', 'Entry': 'acc'})
    phens_df = pd.concat([cos_proteins, adhd, scz, t2d])
    return phens_df


if __name__ == '__main__':
    hc = pd.read_csv(cfg.data['hc'] + '/genetrek-2022-08-23.tsv', sep='\t')
    gene_protein = pd.read_excel(cfg.data['hc'] + '/HC_NDD_genes list_1737to1721.xlsx')
    # merge with protein names
    hc = pd.merge(gene_protein, hc, left_on='Gene_name', right_on='Gene', how='inner')
    hc = hc[['Gene_name', 'acc', 'SysID primary genes (updated 2021-11-18)', 'Epilepsy ALL genesv2',
             'SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)']].drop_duplicates()
    ##
    hc_phens = hc_all_phens_generator(hc)
    hc = pd.merge(hc, hc_phens, on='acc')
    hc = hc.drop(columns=['SysID primary genes (updated 2021-11-18)', 'Epilepsy ALL genesv2',
                          'SFARI_1 Genes (updated 2021 Q3)', 'SFARI_23S Genes (updated 2021 Q3)', 'all_phens'])
    ## concatenating hc and cosminc
    ## cosmic genes downloaded August 28th
    hc = pd.concat([hc, phen_df_maker()])
    hc.to_csv(cfg.data['hc'] + '/smaller-hc-with-phens-column')
