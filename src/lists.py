import pandas as pd
import config as cfg

hc = pd.read_csv(cfg.data['hc'] + '/smaller-hc-with-phens-column')
phens_lst = ['ID', 'ASD', 'ADHD', 'Ep']
hc = hc.loc[hc.phenotype.isin(phens_lst)]
hc = hc.drop(columns=['Unnamed: 0', 'phenotype'])
mobidb = pd.read_csv(cfg.data['clin'] + '/mobidb_result_2022-10-03T14_28_50.889Z.tsv', sep='\t')
features = ['curated-disorder-uniprot', 'curated-disorder-disprot', 'curated-disorder-priority',
            'curated-disorder-ideal', 'curated-disorder-merge']
mobidb_disprot = mobidb.loc[mobidb['feature'] == 'curated-disorder-disprot']
mobidb_disprot = pd.merge(hc, mobidb_disprot, on='acc')

mobidb_ideal = mobidb.loc[mobidb['feature'] == 'curated-disorder-ideal']
mobidb_ideal = pd.merge(hc, mobidb_ideal, on='acc')

mobidb_af = mobidb.loc[mobidb['feature'] == 'prediction-plddt-alphafold']
mobidb_af = pd.merge(hc, mobidb_af, on='acc')
mobidb_af.to_csv(cfg.data['hc'] + '/hc_list_mobidb_alphafold.csv')

disprot_ideal = pd.concat([mobidb_disprot, mobidb_ideal]).sort_values('feature').drop_duplicates(subset=['acc'],
                                                                                                 keep='first')  # 322 unique proteins
phens_lst = ['ID', 'ASD', 'ADHD', 'Ep']
disprot_ideal.to_csv(cfg.data['hc'] + '/hc_list_mobidb_based_on_Disprot_and_IDEAL.csv')

disprot_ideal_af = pd.concat([disprot_ideal, mobidb_af]).sort_values('feature')
disprot_ideal_af.to_csv(cfg.data['hc'] + '/hc_list_mobidb_based_on_Disprot_IDEAL_Alphafold.csv')

pdb = pd.read_csv(cfg.data['pdb'] + '/uniprot-to-pdb.txt', on_bad_lines='skip', sep='  ').reset_index()
pdb = pdb.drop(columns=['level_2', 'Method', 'Unnamed: 2', 'Resolution'])
pdb.columns = ['pdb', 'method', 'resolution', 'acc']
pdb['acc'] = pdb['acc'].str.replace('\(', "")
pdb['acc'] = pdb['acc'].str.replace('\)', "")

all_pdb = pd.merge(disprot_ideal_af, pdb, on='acc')
all_pdb.to_csv(cfg.data['pdb'] + '/Disprot-ideal-alphafold_with_PDB.csv')

## preparing NDD protein list for curation based on List from federica with PDBs added here 11/10/2022
hc_pdb = pd.merge(hc, pdb, on='acc')
curate = pd.read_csv(cfg.data['dis-cu'] + '/Homo_sapiens_9606 - Homo_sapiens_9606.tsv', sep='\t')
curate_hc_pdb = pd.merge(hc_pdb, curate, on='acc')
curate_hc_pdb.to_csv(cfg.data['dis-cu'] + '/hc-proteins-for-curation.csv')
