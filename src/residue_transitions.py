import pandas as pd
import config as cfg
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def residue_heatmapper(df_lst, hmap_title_lst, filename):
    # input df is disorder_majority or filtered_dis_maj if preferred
    new_pivot_df_lst, aa_symbols_lst = heatmap_pivotdf_maker(df_lst)
    print('1')
    sns.set(font_scale=2.4)
    fig, axes = plt.subplots(len(new_pivot_df_lst), 1, figsize=(30, 20 * len(new_pivot_df_lst)))
    print('2')
    for i, (ax, d, t) in enumerate(zip(axes.reshape(-1), new_pivot_df_lst, hmap_title_lst)):
        print('3')
        sb = sns.heatmap(d, cmap="viridis", annot=True, fmt='g', linewidth=1.1, linecolor='black', ax=ax, square=True,
                         cbar_kws={'label': 'Transition percentage'}
                         )
        print('4')
        ax.set_title(t, fontsize=40)
        ax.set_xlabel('Variant residues', fontsize=30)
        ax.set_ylabel('Original residues', fontsize=30)
        ax.tick_params(axis='x', colors='red')
        ax.set_xticklabels(aa_symbols_lst)
        print('5')

        # if i < (len(new_pivot_df_lst) - 1):
        #     sb.set(xticklabels=[])
        #     sb.set(xlabel=None)
    plt.tight_layout()
    plt.savefig(cfg.plots['var-hms'] + '/' + filename + '.png', dpi=120)
    plt.show()

    # ax.set_title(hmap_title, fontsize=40)
    # plt.xlabel('Variant residues', fontsize=25)
    # # ax.xaxis.label.set_color('red')
    # plt.ylabel('Original residues', fontsize=25)
    # # ax.yaxis.label.set_color('blue')
    return


def heatmap_pivotdf_maker(df_lst):
    # todo check if the transitions are counted correctly, because the numbers here are different
    new_pivot_df_lst = []
    for df in df_lst:
        res_df = df[['acc', 'var_id', 'ref_aa', 'alt_aa']]
        aa_categ_order_x = ['Cys', 'Met', 'Val', 'Leu', 'Ile', 'Trp', 'Tyr', 'Phe', 'His', 'Thr', 'Asn', 'Gln', 'Lys',
                            'Arg', 'Asp', 'Glu', 'Ala', 'Ser', 'Gly', 'Pro']
        res_df = res_df.loc[res_df.alt_aa.isin(aa_categ_order_x)]
        residue_ser = res_df.groupby(['ref_aa', 'alt_aa']).size().to_frame(name='size').reset_index()
        residue_pivot_df = pd.pivot_table(residue_ser, values='size', index='ref_aa', columns='alt_aa')
        aa_categ_order_y = aa_categ_order_x.__reversed__()
        residue_pivot_df = residue_pivot_df.reindex(columns=aa_categ_order_x)
        residue_pivot_df = residue_pivot_df.reindex(aa_categ_order_y)
        # max_value = residue_pivot_df.max() # this finds max of each column and the result percentage is comparative
        # only for each column
        max_value = residue_ser['size'].max()
        residue_pivot_df = (residue_pivot_df.div(max_value)).round(3) * 100
        residue_pivot_df = residue_pivot_df.apply(pd.to_numeric)
        new_pivot_df_lst.append(residue_pivot_df)
    return residue_pivot_df, aa_categ_order_x


if __name__ == '__main__':
    mobidb_checked_vars = pd.read_csv(cfg.data['clin'] + '/vars-in_and_out_idr-checked-by-mobidb.csv')
    ## chenge aa 3 letter symbols to one letter

    ##
    all_snv = mobidb_checked_vars.loc[mobidb_checked_vars['Type'] == 'single nucleotide variant']
    in_idr_vars = all_snv.loc[all_snv['isin_idr'] == 1]

    a,b = heatmap_pivotdf_maker([in_idr_vars])
    # sns.heatmap(a)
    fig, ax = plt.subplots(figsize=(15, 15))
    sb = sns.heatmap(a, cmap="plasma", annot=True, fmt='g', linewidth=0.8, linecolor='white', square=True,
                cbar_kws={'label': 'Transition percentage'}, ax=ax, annot_kws={"fontsize":15})
    sb.set_xticklabels(sb.get_xmajorticklabels(), fontsize=18, rotation=45)
    sb.set_yticklabels(sb.get_ymajorticklabels(), fontsize=18, rotation=45)
    sb.tick_params(axis='x', colors='red')
    sb.set_title('SNV residue variations', fontsize=40)
    sb.set_xlabel('Alternate residues', fontsize=25)
    sb.set_ylabel('Reference residues', fontsize=25)
    plt.savefig(cfg.plots['hm'] + '/residue-transition-hmap.png', dpi=120)
    plt.show()

    # residue_heatmapper([in_idr_vars, all_snv], ['SNV-associated Residue transitions - in IDRs',
    #                                             'SNV-associated Residue transitions - All',
    #                                             'Difference (inIDR - Total)'], 'heatmap')
