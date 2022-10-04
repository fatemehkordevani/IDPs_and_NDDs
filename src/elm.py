import pandas as pd
import config as cfg
import matplotlib.pyplot as plt


def df_preparare_multiple_donuts(df):
    df = df[['phenotype', 'ELMType']]
    df_grouped = df.groupby(['phenotype', 'ELMType']).size().reset_index()
    table = pd.pivot_table(df_grouped, index=['phenotype'], columns=['ELMType'])
    # Delete columns with Nan values and convert to integers
    table = table.fillna(0).astype(int)
    # delete numbers from column names
    table.columns = table.columns.get_level_values(1)
    ## changing place of SCZ with T2d for the plot to have last axes with complete legend
    table = table.reset_index()
    b, c = table.iloc[5].copy(), table.iloc[6].copy()
    table.iloc[5], table.iloc[6] = c, b
    table.index = table.phenotype
    del table['phenotype']
    return table


def multiple_donut_plots(df):
    # from: https://sharkcoder.com/data-visualization/mpl-pie-charts
    font_color = '#525252'
    colors = ['#f7ecb0', '#ffb3e6', '#99ff99', '#66b3ff', '#c7b3fb', '#ff6666', '#f9c3b7']

    fig, axes = plt.subplots(3, 3, figsize=(10, 10), facecolor='#e8f4f0')
    fig.delaxes(ax=axes[2, 2])

    for i, (idx, row) in enumerate(df.head(8).iterrows()):
        ax = axes[i // 3, i % 3]
        row = row[row.gt(row.sum() * .01)]
        ax.pie(row,
               labels=row.values,
               startangle=30,
               wedgeprops=dict(width=.5),  # For donuts
               colors=colors,
               textprops={'color': font_color})
        ax.set_title(idx, fontsize=16, color=font_color)

    ax.legend(df.columns, loc='upper right', bbox_to_anchor=(3.5, 0.9), prop={'size': 13})

    fig.subplots_adjust(wspace=.2)  # Space between charts
    fig.delaxes(axes[2][1])

    title = fig.suptitle('ELM classes with variants occurring inside them', y=.95, fontsize=20, color=font_color)
    # To prevent the title from being cropped
    plt.subplots_adjust(top=0.85, bottom=0.15)
    plt.savefig(cfg.plots['dc'] + '/multiple_donuts_ELM_classes_per_phenotype_var_in_elm.png')
    plt.show()
    return


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
    # print(set(transformed_regions))
    return set(transformed_regions)


def isin_elm_array_maker(input_df, checked_col_name):
    # input df is merged df of mobidb and mutation positions from mutinfo
    ## checks if mutation position is in startend disorder region of mobidb or not
    array_is_in = []  # will be filled with boolean of 0,1
    for index, row in input_df.iterrows():
        set_disorder_region = expand_regions(row.startend)  # temp set of data, convert each startend lst to a set,

        if int(row.pr_pos) in set_disorder_region:
            array_is_in.append('1')
        else:
            array_is_in.append('0')
    input_df[checked_col_name] = array_is_in
    return input_df


def var_countcol_creator(df):  # a percentage column for variants in IDR / all Vars
    # input it the out put of df_feature_filtered method
    checked_col_name = 'isin_elm'
    df[checked_col_name] = df[checked_col_name].apply(pd.to_numeric)
    var_count_dic = (df.pivot_table(columns=['acc'], aggfunc='size')).to_dict()
    # dict values to list so that I could append other data to dict values
    for key, value in var_count_dic.items():
        var_count_dic[key] = [value]
    for k in var_count_dic:  # adds the count of each acc var being in idr
        each_acc_vars_in_count = df.loc[df['acc'] == k, checked_col_name].sum()
        var_count_dic[k].append(each_acc_vars_in_count)
    # dict to df
    var_count_df = pd.DataFrame.from_dict(var_count_dic, orient='index', columns=['total_vars', 'in_elm_vars'])
    var_count_df = var_count_df.reset_index()
    var_count_df['out_elm_vars'] = (var_count_df['total_vars'] - var_count_df['in_elm_vars'])
    var_count_df = var_count_df.reset_index()
    var_count_df = var_count_df.rename(columns={'index': 'acc'})
    mrg_var_and_count_df = pd.merge(df, var_count_df, on='acc')
    # mrg_var_and_count_df = mrg_var_and_count_df.drop(columns=['Unnamed: 0', 'Unnamed: 0_x', 'Unnamed: 0_y'])
    return mrg_var_and_count_df


if __name__ == '__main__':

    # downloaded 3rd of october
    # from http://elm.eu.org/instances/?q=*&instance_logic=&taxon=Homo+sapiens&submit=submit&reset_form=Reset
    elm = pd.read_csv(cfg.data['elm'] + '/elm_instances (1).tsv', sep='\t')
    elm = elm.rename(columns={'Primary_Acc': 'acc'})
    hc = pd.read_csv(cfg.data['hc'] + '/smaller-hc-with-phens-column')
    hc_elm = pd.merge(hc, elm, on='acc')
    # multiple_donut_plots(df_preparare_multiple_donuts(hc_elm))
    ## clinvar df with var positions related to our phenotypes
    elm['startend'] = elm['Start'].astype('string') + '..' + elm['End'].astype('string')
    elm = elm.drop(columns=['ProteinName', 'Accessions', 'Start', 'End', 'Methods', 'InstanceLogic', 'Organism'])
    elm = elm.groupby(['ELMType', 'acc'])['startend'].apply(','.join).reset_index()

    hc_clin = pd.read_csv(cfg.data['hc'] + '/clinvar_and_phenotypes_mutation_positions.csv')
    clin_elm = pd.merge(hc_clin, elm, on='acc')
    clin_elm = var_countcol_creator(isin_elm_array_maker(clin_elm, 'isin_elm'))
    in_elm_vars = clin_elm.loc[clin_elm['isin_elm'] == 1]

    multiple_donut_plots(df_preparare_multiple_donuts(in_elm_vars))
