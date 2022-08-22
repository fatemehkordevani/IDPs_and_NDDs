from io import StringIO
import pandas as pd
import config as cfg


# ndd_clinvar = pd.read_csv(cfg.data['clin'] + '/ndd-clin-pos-merged.tsv') annovar = pd.read_csv(cfg.data['clin'] +
# '/query.output.exome_summary.csv') annovar['gnomAD_exome_ALL'] = annovar.loc[(annovar['gnomAD_exome_ALL'] != '.') &
# (annovar['gnomAD_exome_ALL'] != 'nan'), 'gnomAD_exome_ALL'].astype(float).apply(lambda x: '%.15f' % x)
# annovar.to_csv(cfg.data['clin'] + '/annovar-ndd-vars.tsv')

## from here started at Monday Aug 8th
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


def mutidr_bool_array_maker(input_df):
    # input df is merged df of mobidb and mutation positions from mutinfo
    ## checks if mutation position is in startend disorder region of mobidb or not
    array_is_in = []  # will be filled with boolean of 0,1
    for index, row in input_df.iterrows():
        set_disorder_region = expand_regions(row.startend)  # temp set of data, convert each startend lst to a set,

        if int(row.pr_pos) in set_disorder_region:
            array_is_in.append('1')
        else:
            array_is_in.append('0')
    input_df['isin_idr'] = array_is_in
    return input_df


# ## clinvar version: variant summary.txt 135,193,289	2022-08-01 17:10:07	2022-08-01 17:10:07
# clinvar = pd.read_csv(cfg.data['clin'] + '/variant_summary.txt', sep='\t', low_memory=False)
# clinvar = clinvar.drop(columns=['Assembly', 'ChromosomeAccession', 'Start', 'Stop', 'PositionVCF', 'ReferenceAlleleVCF',
#                                 'AlternateAlleleVCF'])
# clinvar = clinvar.drop_duplicates()
# clinvar = clinvar.loc[(clinvar['ClinicalSignificance'] == 'Pathogenic') |
#                       (clinvar['ClinicalSignificance'] == 'Likely pathogenic')
#                       | (clinvar['ClinicalSignificance'] == 'Pathogenic/Likely pathogenic')]
# ## Gene4denovo candidate genes   Release version 1.2 (Updated: 07/08/2022)
# candidates = pd.read_excel(cfg.data['gene4'] + '/Candidate_gene_1.2.xlsx', engine='openpyxl', usecols=['Groups', 'Gene_symbol', 'FDR'])
# candidates = candidates.loc[candidates['FDR'] < 0.05]
# cand_genes = candidates['Gene_symbol'].unique().tolist()
# ## the protein IDs corresponding to candidate genes are retrieved from uniprot on 09/08/2022
# ## you should figure out genes with more than one protein
# cand_genes_pr_ids = pd.read_csv(cfg.data['clin'] + '/uniprot-cand-genes_protein_IDs.tsv', sep='\t')
# del cand_genes_pr_ids['Reviewed']
# ## merging to also have g4dn phenotype names
# cand_genes_pr_ids = pd.merge(candidates, cand_genes_pr_ids, left_on='Gene_symbol', right_on='Gene_name', how='inner')
# ## only clinvar row that are representing our candidate genes, ## merging protein IDs with clinvar
# clinvar_ndd = clinvar.loc[clinvar.GeneSymbol.isin(cand_genes)]
# clinvar_ndd = pd.merge(cand_genes_pr_ids, clinvar_ndd, left_on='Gene_name', right_on='GeneSymbol', how='inner')
# clinvar_ndd.to_csv(cfg.data['clin'] + '/clinvar-ndd-candidate-genes.csv')
clinvar_ndd = pd.read_csv(cfg.data['clin'] + '/clinvar-ndd-candidate-genes.csv')
del clinvar_ndd['Unnamed: 0']
del clinvar_ndd['Gene_name']
##
clinvar_ndd[['Name', 'pr_change']] = clinvar_ndd['Name'].str.split('p.', 1, expand=True)
# clinvar_ndd[['pr_change']] = clinvar_ndd['pr_change'].str.split('p.', 1, expand=True)[1]
clinvar_ndd = clinvar_ndd.dropna(subset=['pr_change'])
clinvar_ndd['pr_change'] = clinvar_ndd['pr_change'].str.replace(r')', '')
clinvar_ndd = clinvar_ndd.loc[~clinvar_ndd.pr_change.str.contains('=')]
clinvar_ndd['ref_aa'] = clinvar_ndd['pr_change'].str[:3]
clinvar_ndd['alt_aa'] = clinvar_ndd['pr_change'].str[-3:]
clinvar_ndd['pr_pos'] = clinvar_ndd.pr_change.str.extract('(\d+)')
clinvar_ndd = clinvar_ndd.dropna(subset=['pr_pos'])
## change three letter aa symbol to one letter
clinvar_ndd = clinvar_ndd.loc[clinvar_ndd['PhenotypeList'] != 'not provided']
phens_lst = clinvar_ndd['PhenotypeList'].unique().tolist()
with open(cfg.data['clin'] + '/your_file.txt', 'w') as f:
    for i in phens_lst:
        f.write("%s\n" % i)
## keep only phenotypes with my keywords
## what about capital or small letters?
# phens_keywords = ['Intellectual', 'development', 'motor', 'mobility', 'movement', 'autism', 'spectrum',
#                   'attention deficit', 'seizure', 'retardation', 'syndrome', 'brain', 'epileptic', 'epilepsy', 'schizo']
phens_keywords2 = ['intellectual disability', 'developmental delay', 'neurodevelopmental disorder', 'motor', 'mobility',
                   'movement', 'autis', 'attention deficit', 'retard', 'epilep', 'schizo']
phen_keyword_dict = {'ID': ['Intellectual disability', 'INTELLECTUAL DEFICIENCY'],
                     'DD': ['developmental delay', 'neurodevelopmental delay'], 'NDD': ['neurodevelopmental disorder'],
                     'Motor disorders': ['motor, movement ,mobility'], 'ASD': ['autism spectrum disorder'],
                     'ADHD': ['attention deficit'],
                     'EE': ['epilepsy', 'epileptic encephalopathy'], 'SCZ': ['schizophrenia']}
phens_filtered_list = []
for kword in phens_keywords2:
    for phen in phens_lst:
        if kword.upper() in phen.upper():
            phens_filtered_list.append(phen)
        else:
            continue
## replacing phens_filtered_list items with phen_keyword_dict, and then exploding the df
phens_short = []
p = 0
for phen_orig in phens_filtered_list:
    phen_orig = phen_orig.upper()
    for k, v in phen_keyword_dict.items():
        for each_v in v:
            each_v = each_v.upper()
            if each_v in phen_orig:
                phen_orig = phen_orig.replace(each_v, k)
    print(phen_orig)
    phens_short.append(phen_orig)
    p += 1
    print(len(phens_short))
phens_kwd_tuple = list(zip(phens_filtered_list, phens_short))
## filter clinvar_ndd for only my phenotypes
for phen_first in phens_filtered_list:
    clinvar_ndd.loc[clinvar_ndd['PhenotypeList'] == phen_first, 'phen'] = phens_short[phens_filtered_list.index(phen_first)]
clinvar_ndd = clinvar_ndd[clinvar_ndd['phen'].notna()]
# ## Mobidb downloaded on Aug 8th 2022, 78,106 entries (but 77629 prs after filtering for disorder consensus)
mobidb = pd.read_csv(cfg.data['clin'] + '/mobidb_result_2022-08-08T13_12_39.379Z.tsv', sep='\t')
mobidb = mobidb.rename(columns={'start..end': 'startend'})
mobidb = mobidb.loc[mobidb['feature'] == 'prediction-disorder-th_50']  # (77629, 6)
##
clin_mobidb_ndd = pd.merge(clinvar_ndd, mobidb, on='acc', how='inner')  # (35260, 37)
## check if modification position is in idr
clin_mobidb_ndd = mutidr_bool_array_maker(clin_mobidb_ndd)
phenotypes = ['ASD', 'EE', 'ID', 'DD', 'SCZ', 'NDDs', 'Control', 'Complex motor steriotypies']
# # should we extract control proteins from other phenotypes? even though they do not exist here independantly?
# that is weird by the way!
clin_mobidb_ndd = clin_mobidb_ndd.loc[clin_mobidb_ndd.Groups.isin(phenotypes)]
in_idr_vars = clin_mobidb_ndd.loc[clin_mobidb_ndd['isin_idr'] == '1']
