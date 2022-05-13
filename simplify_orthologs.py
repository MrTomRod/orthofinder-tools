import os
from orthofinder_tools.orthogroup_to_gene_name import OrthogroupToGeneName


def simplify_orthologs(orthofinder_tsv: str, fasta_dir: str, out_dir: str, n0: bool = True, file_endings='fasta'):
    assert os.path.isdir(out_dir), f'Does not exist: {out_dir}'

    otg = OrthogroupToGeneName(
        fasta_dir=fasta_dir,
        file_endings=file_endings
    )

    if n0:
        otg.load_hog(
            n0_tsv=orthofinder_tsv
        )
    else:
        otg.load_og(
            og_tsv=orthofinder_tsv
        )

    prefix = 'HOG' if n0 else 'OG'

    otg.save_orthogroup_to_gene_ids(F'{out_dir}/{prefix}_gene_ids.tsv')
    otg.save_majority_df(F'{out_dir}/{prefix}_majority_df.tsv')
    otg.save_orthogroup_to_best_name(F'{out_dir}/{prefix}_best_name.tsv')


if __name__ == "__main__":
    import fire  # pip install fire # automated argparse

    fire.Fire(simplify_orthologs)
