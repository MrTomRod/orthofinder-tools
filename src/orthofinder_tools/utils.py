import os
import pandas as pd

ALLOWED_RESULT_TYPES = ['gene-list', 'count', 'boolean']


def get_apply_fn(result_type: str):
    if result_type == 'gene-list':
        apply_fn = lambda x: [] if pd.isnull(x) else x.split(', ')
    elif result_type == 'count':
        apply_fn = lambda x: 0 if pd.isnull(x) else x.count(', ') + 1
    elif result_type == 'boolean':
        apply_fn = lambda x: not pd.isnull(x)
    else:
        raise AssertionError(f'Failed to understand result_type={result_type}. Must be one of {ALLOWED_RESULT_TYPES}.')
    return apply_fn


def load_og(og_tsv: str, result_type: str = 'gene-list'):
    og_tsv = os.path.abspath(os.path.expanduser(og_tsv))
    assert os.path.isfile(og_tsv), F'og_tsv does not exist: "{og_tsv}"'

    apply_fn = get_apply_fn(result_type)

    # read "Orthogroups.tsv"
    gene_ids_df = pd.read_csv(og_tsv, sep='\t', dtype=str)

    # sanity checks
    assert 'Orthogroup' in gene_ids_df.columns, \
        f'The file {og_tsv} does not appear to be a valid Orthogroup.tsv: ' \
        f'The column "Orthogroup" is missing. Try setting --hog=False'

    gene_ids_df.set_index('Orthogroup', inplace=True)

    gene_ids_df = gene_ids_df.map(apply_fn)

    return gene_ids_df


def load_hog(hog_tsv: str, result_type: str = 'gene-list'):
    hog_tsv = os.path.abspath(os.path.expanduser(hog_tsv))
    assert os.path.isfile(hog_tsv), F'hog_tsv does not exist: "{hog_tsv}"'

    apply_fn = get_apply_fn(result_type)

    # read "N0.tsv" and drop extra columns
    gene_ids_df = pd.read_csv(hog_tsv, sep='\t', dtype=str)

    # sanity checks
    for col in ['HOG', 'OG', 'Gene Tree Parent Clade']:
        assert col in gene_ids_df.columns, f'The file {hog_tsv} does not appear to be a valid N0.tsv: ' \
                                           f'The column "{col}" is missing. Try setting --hog=False'

    gene_ids_df.set_index('HOG', inplace=True)

    gene_ids_df.drop(columns=['OG', 'Gene Tree Parent Clade'], inplace=True)

    gene_ids_df = gene_ids_df.map(apply_fn)

    return gene_ids_df
