import pandas as pd
import config as cfg
import matplotlib.pyplot as plt
import seaborn as sns


def var_llps_df_merger(source):
    # df with checked vars for being in idr or not (based on disorder majority consensus)
    # disorder_majority = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv',
    #                                 usecols=
    #                                 ['acc', 'var_id', 'orig_aa', 'var_aa', 'position', 'isin_idr', 'total_vars',
    #                                  'in_idr_vars', 'out_idr_vars'])
    # in_idr_var_df = disorder_majority.loc[disorder_majority['isin_idr'] == 1]
    # in_idr_var_lst = in_idr_var_df['acc'].unique().tolist()
    # NDD (in IDR vaqriants df but just for NDD protein)
    clin_mobidb_ndd = pd.read_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb.csv')
    ndd_invar_df = clin_mobidb_ndd.loc[clin_mobidb_ndd['isin_idr'] == 1]
    # LLPS
    phasepro_df, phasepro_lst = phase_pro()
    phasep_df, phasep_lst = phasep()
    # merging
    if source == 'phasepro':
        # var_phasepro_mrg = pd.merge(in_idr_var_df, phasepro_df, on='acc')
        ndd_var_phasepro_mrg = pd.merge(ndd_invar_df, phasepro_df, on='acc')
        return  ndd_var_phasepro_mrg
    elif source == 'phasep':
        # var_phasep_mrg = pd.merge(in_idr_var_df, phasep_df, left_on='acc', right_on='UniprotEntry')
        ndd_var_phasep_mrg = pd.merge(ndd_invar_df, phasep_df, left_on='acc', right_on='UniprotEntry')
        return ndd_var_phasep_mrg


def phase_pro():
    phasepro = pd.read_csv(cfg.data['fp'] + '/phasepro.tsv', sep='\t')
    phasepro.columns = ['common_name', 'name', 'acc', 'organism', 'sequence', 'gene', 'taxon', 'id', 'segment',
                        'boundaries', 'region', 'partners', 'determinants', 'forms', 'organelles', 'pmids',
                        'description',
                        'experiment_llps', 'in_vivo', 'in_vitro', 'experiment_state', 'rna_req', 'ptm_affect',
                        'disease',
                        'splice', 'interaction', 'membrane_cluster', 'partner_dep', 'rna_dep', 'ptm_dep',
                        'domain-motif_interaction', 'discrete_oligo', 'under_annote', 'annotator', 'functional_class',
                        'date']
    phasepro_lst = phasepro['acc'].unique().tolist()
    return phasepro, phasepro_lst


def phasep():
    phasep = pd.read_excel(cfg.data['fs'] + '/High throughput Data V1.3.xlsx', engine='openpyxl')
    phasep = phasep.drop(['No', 'protein material states', 'Mutation/disease', 'Organism'], axis=1)
    phasep_lst = phasep['UniprotEntry'].unique().tolist()
    return phasep, phasep_lst


if __name__ == '__main__':
    ## the two llps datasets have 28 proteins in common
    ndd_subdf = pd.read_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb.csv')
    ndd_prs = ndd_subdf.loc[(ndd_subdf['phenotype'] == 'ASD') | (ndd_subdf['phenotype'] == 'ID') | (ndd_subdf['phenotype'] == 'Ep')]
    ndd_prs = ndd_prs['acc'].unique().tolist()


    ndd_varin_phasep_df = var_llps_df_merger('phasep')  # 2296 # 840
    ndd_varin_phasepro_df = var_llps_df_merger('phasepro')  # 134 # 30
    # lst of NDD proteins with var in IDR and involved in LLPS
    ndd_varin_phasep_lst = ndd_varin_phasep_df['acc'].unique().tolist()  # 41
    ndd_varin_phasepro_lst = ndd_varin_phasepro_df['acc'].unique().tolist()  # 5
    ## new mlo dis based on my ndds, var in idr ndds, merged two llps dfs
    ndd_llps_merged_lst = ndd_varin_phasep_lst + ndd_varin_phasepro_lst

    ## Disphase
    disphase_pr_lst = pd.read_csv(cfg.data['fs-df'] + '/disphase_pr-lst.csv')
    disphase_pr_lst = disphase_pr_lst['uniprot_acc'].unique().tolist()
    ## NDD intersect with disphase
    ndd_disphase = list(set(ndd_prs) & set(disphase_pr_lst))
    # print(','.join(ndd_disphase))

    ## DrLLPS
    drllps = pd.read_excel(cfg.data['fs-dr'] + '/Table S1.xlsx')
    drllps.columns = drllps.iloc[0]
    drllps_prs = drllps['UniProt ID'].unique().tolist()
    disph_drllps = list(set(drllps_prs) & set(ndd_disphase))
    drllps = drllps.loc[drllps['UniProt ID'].isin(disph_drllps)]
    drllps_types = drllps.groupby('LLPS Type').count()


    ## file from this paper:https://www.sciencedirect.com/science/article/pii/S2001037021002804
    ## challenge, open this file correctly, for now not very necessary
    # mlo = pd.read_excel(cfg.data['mlo']+'/complementary-data-paper-S2001037021002804.xlsx', engine='openpyxl')
    ## STRING analysis of ndd_llps_merged_lst with 46 total proteins
    # molecular_function = pd.read_csv(cfg.data['fs-str'] + '/enrichment.Function.tsv', sep='\t')
    # cell_process = pd.read_csv(cfg.data['fs-str'] + '/enrichment.Process.tsv', sep='\t')
    # count number of variants in llps dataset
    llps_ndd_dismaj = ndd_subdf.loc[ndd_subdf.acc.isin(ndd_disphase)]
    llps_ndd_dismaj = llps_ndd_dismaj.drop_duplicates()
    llps_ndd_dismaj = llps_ndd_dismaj.rename(columns={"in_idr_vars_perc": "IDR mutation fraction (NDD-associated Proteins with LLPS roles)"})
    # g = sns.violinplot(x=llps_ndd_dismaj['IDR mutation fraction (NDD-associated Proteins with LLPS roles)'])
    # g.set_xlim(0, 1)
    # plt.savefig(cfg.plots['vp'] + '/ndd-llps-roles-disphase.png')
    #
    # plt.show()
    #
    # ndd_subdf = ndd_subdf.rename(columns={"in_idr_vars_perc": "IDR mutation fraction (All NDD-associated proteins)"})
    # g = sns.violinplot(x=ndd_subdf['IDR mutation fraction (All NDD-associated proteins)'])
    # g.set_xlim(0, 1)
    # plt.savefig(cfg.plots['vp'] + '/all-ndd-HC.png')
    # plt.show()
    mlo = pd.read_csv(cfg.data['mlo-d'] + '/mlodisdb_components.csv')
    mlo = mlo.loc[mlo.Entry.isin(ndd_llps_merged_lst)]
    mlo_count = mlo.groupby('MLO').count()
    mlo_count = mlo_count.reset_index()
    mlo_count = mlo_count.drop(columns=['Entry name', 'Gene names', 'Source'])
    mlo_count = mlo_count.rename(columns={'Entry': 'Count'})



