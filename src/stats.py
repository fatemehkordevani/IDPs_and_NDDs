import pandas as pd
import numpy as np
import config as cfg
import scipy.stats as stats
from pprint import pprint

if __name__ == '__main__':
    ## Stat test for disorder content based on mobidb
    disorder_stat = pd.read_csv(cfg.data['stats'] + '/mobidb-disorder-for-stats')
    disorder_stat = disorder_stat.fillna(0)
    fvalue, pvalue = stats.f_oneway(disorder_stat['Human'], disorder_stat['Brain'], disorder_stat['ASD'],
                                    disorder_stat['Ep'], disorder_stat['ID'], disorder_stat['ADHD'],
                                    disorder_stat['SCZ'], disorder_stat['Cancer'], disorder_stat['T2D'])
    print(fvalue, pvalue)
    ## Stat test for LIP content based on mobidb
    lip_stat = pd.read_csv(cfg.data['stats'] + '/mobidb-lip-for-stats')
    lip_stat = lip_stat.fillna(0)
    fvalue, pvalue = stats.f_oneway(lip_stat['Human'], lip_stat['Brain'], lip_stat['ASD'],
                                    lip_stat['Ep'], lip_stat['ID'], lip_stat['ADHD'],
                                    lip_stat['SCZ'], lip_stat['Cancer'], lip_stat['T2D'])
    print(fvalue, pvalue)
    ## Stat test for Domain content based on mobidb
    domain_stat = pd.read_csv(cfg.data['stats'] + '/mobidb-domain-for-stats')
    domain_stat = domain_stat.fillna(0)
    fvalue, pvalue = stats.f_oneway(domain_stat['Human'], domain_stat['Brain'], domain_stat['ASD'],
                                    domain_stat['Ep'], domain_stat['ID'], domain_stat['ADHD'],
                                    domain_stat['SCZ'], domain_stat['Cancer'], domain_stat['T2D'])
    print(fvalue, pvalue)
    # Stat for disorder content based on disprot
    dis_dis = pd.read_csv(cfg.data['dis'] + '/disprot_disorder_for_dis_content_boxplot_r.csv')
    dis_dis = dis_dis.pivot(index=['index', 'acc'], columns='phenotype', values='content_fraction')
    dis_dis = dis_dis.fillna(0)
    fvalue, pvalue = stats.f_oneway(dis_dis['ASD'],dis_dis['Ep'], dis_dis['ID'], dis_dis['ADHD'],
                                    dis_dis['SCZ'], dis_dis['Cancer'], dis_dis['T2D'])
    print(fvalue, pvalue)
    # Stat for protein binding based on Disprot
    dis_bind = pd.read_csv(cfg.data['dis'] + '/disprot_binding_for_dis_content_boxplot_r.csv')
    dis_bind = dis_bind.pivot(index=['index', 'acc'], columns='phenotype', values='content_fraction')
    dis_bind = dis_bind.fillna(0)
    fvalue, pvalue = stats.f_oneway(dis_bind['ASD'],dis_bind['Ep'], dis_bind['ID'], dis_bind['ADHD'],
                                    dis_bind['SCZ'], dis_bind['Cancer'], dis_bind['T2D'])
    print(fvalue, pvalue)
    # Stat for disorder to order transition based on Disprot
    dis_ord = pd.read_csv(cfg.data['dis'] + '/disprot_D-to-O_for_dis_content_boxplot_r.csv')
    dis_ord = dis_ord.pivot(index=['index', 'acc'], columns='phenotype', values='content_fraction')
    dis_ord = dis_ord.fillna(0)
    fvalue, pvalue = stats.f_oneway(dis_ord['ASD'], dis_ord['Ep'], dis_ord['ID'], dis_ord['ADHD'],
                                    dis_ord['SCZ'], dis_ord['Cancer'], dis_ord['T2D'])
    print(fvalue, pvalue)
