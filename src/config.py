import os

absolute = os.path.abspath(os.getcwd())
absolute = absolute + '/../'
data = {
    'gene4': absolute + 'data/gene4denovo',
    'hc': absolute + 'data/hcgenes',
    'rseq': absolute + 'data/refseq',
    '': absolute + 'data',
    'brain': absolute + 'data/brain',
    'csv': absolute + 'data/to-csv',
    'dndb': absolute + 'data/denovodb',
    'phens': absolute + 'data/mobibool-phens',
    'phens-fdr': absolute + 'data/mobibool-phens/acc-phens-FDR',
    'xml': absolute + 'data/xml',
    'xml-p': absolute + 'data/xml/parsed',
    'vars': absolute + 'data/variants',
    'desc-cf': absolute + 'data/description/cf',
    'desc-cc': absolute + 'data/description/cc',
    'desc-len': absolute + 'data/description/length',
    'all-var-desc-cf': absolute + 'data/variants/desc-allvariant/cf',
    'all-var-desc-cc': absolute + 'data/variants/desc-allvariant/cc',
    'all-var-desc-len': absolute + 'data/variants/desc-allvariant/length',
    'dis': absolute + 'data/Disprot',
    'jsn': absolute + 'data/json',
    'ptm': absolute + 'data/ptms',
    'ptm-u': absolute + 'data/ptms/uniprot-ptms',
    'fp': absolute + 'data/llps/phasepro',
    'fs': absolute + 'data/llps/phasep',
    'mlo-d': absolute + 'data/llps/mlo-dis',
    'mlo': absolute + 'data/llps/insight-into-llps',
    'fs-str': absolute + 'data/llps/string',
    'clin': absolute + 'data/clinvar',
    'cos': absolute + 'data/cosmic',
    'go': absolute + 'data/string-go',

}
plots = {
    '': absolute + 'plots',
    'bar': absolute + 'plots/barplot',
    'bar-sptm': absolute + 'plots/barplot/sigpep-transmemb',
    'hm' : absolute + 'plots/heatmap',
    'bp' : absolute + 'plots/boxplot',
    'bp-cf': absolute + 'plots/boxplot/cf',
    'bp-cc': absolute + 'plots/boxplot/cc',
    'bp-len': absolute + 'plots/boxplot/length',
    'vp' : absolute + 'plots/violinplot',
    'vp-cf': absolute + 'plots/violinplot/cf',
    'vp-cc': absolute + 'plots/violinplot/cc',
    'vp-len': absolute + 'plots/violinplot/length',
    'go': absolute + 'plots/go-tables',

}
