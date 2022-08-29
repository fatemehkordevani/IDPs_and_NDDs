import pandas as pd
import config as cfg
import matplotlib.pyplot as plt
import numpy as np


def variant_in_out_idr_stacked_df(df):
    in_idr_df = df.loc[df['isin_idr'] == 1]
    out_idr_df = df.loc[df['isin_idr'] == 0]
    phens_lst = ['ASD', 'ID', 'Ep', 'Cancer']
    phens_count_d = {}
    for phen in phens_lst:
        ## below lddt 50
        in_idr = len(in_idr_df.loc[(in_idr_df['phenotype'] == phen) , 'var_id'])
        out_idr = len(out_idr_df.loc[(out_idr_df['phenotype'] == phen) , 'var_id'])
        phens_count_d[phen] = [in_idr, out_idr]
    df_output = pd.DataFrame(phens_count_d, index= ['Mutations inside IDRs', 'Mutations outside IDRs'])
    df_output = df_output.transpose()
    return df_output


def stacked_percentage_barplotter(df, fig_name):
    # creates stacked barplot with percentages from categorized binding df,
    df_stack = df.apply(lambda x: x * 100 / sum(x), axis=1)
    df_stack.plot(kind='bar', figsize=(25, 15), stacked=True, colormap='RdYlBu').tick_params(axis='both',
                                                                                                labelsize=16)
    # from: towardsdatascience.com/100-stacked-charts-in-python-6ca3e1962d2b
    # and stackoverflow.com/questions/51495982/display-totals-and-percentage-in-stacked-bar-chart-using-dataframe-plot
    for n, x in enumerate([*df.index.values]):
        for percent, count, y_loc in zip(df_stack.loc[x], df.loc[x], df_stack.loc[x].cumsum()):
            # this if chooses what percentages to be shown on bar pieces
            if percent >= 2:
                plt.text(x=n - 0.17,
                         y=(y_loc - percent) + (percent / 2),
                         s=f'{count}\n({np.round(percent, 1)}%)',
                         color="white",
                         fontsize=24,
                         fontweight="bold")
                         # rotation='horizontal')
    plt.legend(bbox_to_anchor=(1.01, 1.02), borderaxespad=0, prop={'size': 16})
    plt.xlabel('Phenotypes', fontsize=26)
    plt.ylabel('Mutation distribution', fontsize=26)
    plt.xticks(fontsize=26, rotation=0)
    plt.yticks(fontsize=26, rotation=0)
    plt.legend(loc=1, prop={'size': 22})
    plt.suptitle('Distribution of mutations within and outside of IDRs', fontsize=30)
    plt.tight_layout()
    plt.savefig(cfg.plots['bar'] + '/' + fig_name + '.png')
    # plt.show()
    return


if __name__ == '__main__':
    ## idr mutations
    clin_mobidb_ndd = pd.read_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb.csv')
    in_idr_vars = clin_mobidb_ndd.loc[clin_mobidb_ndd['isin_idr'] == 1]
    ## lip mutations
    ndd_mobidb_lip = pd.read_csv(cfg.data['hc'] + '/mobidb_vars_in_lips_checked.csv')

    ## stacked bar plot of distribution of variants inside and outside idrs
    stacked_percentage_barplotter(variant_in_out_idr_stacked_df(clin_mobidb_ndd), 'distribution of vars in out idrs')

    # modification_types_freq = in_idr_vars['Type'].value_counts().to_frame().reset_index()
    # modification_types_freq = modification_types_freq.rename(columns={'Type': 'Frequency', 'index': 'Modification types'})
    # modification_types_freq = modification_types_freq.loc[(modification_types_freq['Modification types'] != 'protein only') & (modification_types_freq['Modification types'] != 'Inversion')]
    # plt.bar(modification_types_freq['Modification types'].tolist(), modification_types_freq['Frequency'].tolist())
    # # plt.yscale("log")
    # plt.xticks(modification_types_freq['Modification types'].tolist(), rotation='90', ha='center')
    # plt.title('Frequency of modification types among variants that occur in IDRs')
    # plt.tight_layout()
    # plt.savefig(cfg.plots['bar'] + '/mutation_types1.png')
    # plt.show()

