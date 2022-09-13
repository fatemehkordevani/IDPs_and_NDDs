import pandas as pd
import config as cfg

binding_mo = pd.read_csv(cfg.data['elm'] + '/elm_interaction_domains.tsv', sep = '\t')
