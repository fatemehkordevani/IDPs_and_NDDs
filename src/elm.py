import pandas as pd
import config as cfg

# downloaded 3rd of october
# from http://elm.eu.org/instances/?q=*&instance_logic=&taxon=Homo+sapiens&submit=submit&reset_form=Reset
elm = pd.read_csv(cfg.data['elm'] + '/elm_instances (1).tsv', sep='\t')

