import pandas as pd
import config as cfg
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def variant_in_out_idr_stacked_df(df):
    in_idr_df = df.loc[df['isin_idr'] == 1]
    out_idr_df = df.loc[df['isin_idr'] == 0]
    phens_lst = ['ASD', 'Ep', 'ID', 'ADHD', 'SCZ', 'Cancer', 'T2D']
    phens_count_d = {}
    for phen in phens_lst:
        ## below lddt 50
        in_idr = len(in_idr_df.loc[(in_idr_df['phenotype'] == phen) , 'var_id'])
        out_idr = len(out_idr_df.loc[(out_idr_df['phenotype'] == phen) , 'var_id'])
        phens_count_d[phen] = [in_idr, out_idr]
    df_output = pd.DataFrame(phens_count_d, index= ['Mutations inside IDRs', 'Mutations outside IDRs'])
    df_output = df_output.transpose()
    return df_output


def add_p_value_annotation(fig, array_columns, subplot=None,
                           _format=dict(interline=0.07, text_height=1.07, color='black')):
    ''' Adds notations giving the p-value between two box plot data (t-test two-sided comparison)

    Parameters:
    ----------
    fig: figure
        plotly boxplot figure
    array_columns: np.array
        array of which columns to compare
        e.g.: [[0,1], [1,2]] compares column 0 with 1 and 1 with 2
    subplot: None or int
        specifies if the figures has subplots and what subplot to add the notation to
    _format: dict
        format characteristics for the lines

    Returns:
    -------
    fig: figure
        figure with the added notation
    '''
    # Specify in what y_range to plot for each pair of columns
    y_range = np.zeros([len(array_columns), 2])
    for i in range(len(array_columns)):
        y_range[i] = [1.01 + i * _format['interline'], 1.02 + i * _format['interline']]

    # Get values from figure
    fig_dict = fig.to_dict()

    # Get indices if working with subplots
    if subplot:
        if subplot == 1:
            subplot_str = ''
        else:
            subplot_str = str(subplot)
        indices = []  # Change the box index to the indices of the data for that subplot
        for index, data in enumerate(fig_dict['data']):
            # print(index, data['xaxis'], 'x' + subplot_str)
            if data['xaxis'] == 'x' + subplot_str:
                indices = np.append(indices, index)
        indices = [int(i) for i in indices]
        print((indices))
    else:
        subplot_str = ''

    # Print the p-values
    for index, column_pair in enumerate(array_columns):
        if subplot:
            data_pair = [indices[column_pair[0]], indices[column_pair[1]]]
        else:
            data_pair = column_pair

        # Mare sure it is selecting the data and subplot you want
        # print('0:', fig_dict['data'][data_pair[0]]['name'], fig_dict['data'][data_pair[0]]['xaxis'])
        # print('1:', fig_dict['data'][data_pair[1]]['name'], fig_dict['data'][data_pair[1]]['xaxis'])

        # Get the p-value
        pvalue = stats.ttest_ind(
            fig_dict['data'][data_pair[0]]['y'],
            fig_dict['data'][data_pair[1]]['y'],
            equal_var=False,
        )[1]
        if pvalue >= 0.05:
            symbol = 'ns'
        elif pvalue >= 0.01:
            symbol = '*'
        elif pvalue >= 0.001:
            symbol = '**'
        else:
            symbol = '***'
        # Vertical line
        fig.add_shape(type="line",
                      xref="x" + subplot_str, yref="y" + subplot_str + " domain",
                      x0=column_pair[0], y0=y_range[index][0],
                      x1=column_pair[0], y1=y_range[index][1],
                      line=dict(color=_format['color'], width=2, )
                      )
        # Horizontal line
        fig.add_shape(type="line",
                      xref="x" + subplot_str, yref="y" + subplot_str + " domain",
                      x0=column_pair[0], y0=y_range[index][1],
                      x1=column_pair[1], y1=y_range[index][1],
                      line=dict(color=_format['color'], width=2, )
                      )
        # Vertical line
        fig.add_shape(type="line",
                      xref="x" + subplot_str, yref="y" + subplot_str + " domain",
                      x0=column_pair[1], y0=y_range[index][0],
                      x1=column_pair[1], y1=y_range[index][1],
                      line=dict(color=_format['color'], width=2, )
                      )
        ## add text at the correct x, y coordinates
        ## for bars, there is a direct mapping from the bar number to 0, 1, 2...
        fig.add_annotation(dict(font=dict(color=_format['color'], size=14),
                                x=(column_pair[0] + column_pair[1]) / 2,
                                y=y_range[index][1] * _format['text_height'],
                                showarrow=False,
                                text=symbol,
                                textangle=0,
                                xref="x" + subplot_str,
                                yref="y" + subplot_str + " domain"
                                ))
    return fig


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
    plt.ylabel('Variation distribution', fontsize=26)
    plt.xticks(fontsize=26, rotation=0)
    plt.yticks(fontsize=26, rotation=0)
    plt.legend(loc=1, prop={'size': 22})
    plt.suptitle('Distribution of mutations within and outside of Binding regions', fontsize=30)

    plt.tight_layout()
    plt.savefig(cfg.plots['bar'] + '/' + fig_name + '.png')
    # plt.show()
    return


if __name__ == '__main__':
    ## idr mutations
    clin_mobidb_ndd = pd.read_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb.csv')
    vus_mobidb_ndd = pd.read_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb-VUS.csv')
    ## lip mutations
    ndd_mobidb_lip = pd.read_csv(cfg.data['hc'] + '/mobidb_vars_in_lips_checked.csv')
    vus_mobidb_lip = pd.read_csv(cfg.data['hc'] + '/mobidb_vars_in_lips_checked-VUS.csv')
    ## Disprot disorder variants
    dis_patho = pd.read_csv(cfg.data['dis'] + '/Disprot-disorder-variants-in-out-checked.csv')
    # dis_vus = pd.read_csv(cfg.data['dis'] + '/Disprot-disorder-VUSvariants-in-out-checked.csv')
    ## Disprot binding variants
    dis_bind = pd.read_csv(cfg.data['dis'] + '/Disprot-binding-VUS_variants-in-out-checked.csv')

    ##
    a = variant_in_out_idr_stacked_df(clin_mobidb_ndd)
    # a.to_csv(cfg.data['clin'] + '/stacked_df_patho_var_in_out.csv') # to be imported in Rstudio
    a = a.to_numpy()
    # chi2, p, _, _ = stats.chi2_contingency(a)
    # r, p = stats.pearsonr(a['Mutations inside IDRs'], a['Mutations outside IDRs'])
    # r, p = stats.wilcoxon(a)
    # print(r, p )

    ## stacked bar plot of distribution of variants inside and outside idrs
    # stacked_percentage_barplotter(variant_in_out_idr_stacked_df(dis_bind), 'Disprot-Variants_VUS-binding')

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

