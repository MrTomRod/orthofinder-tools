import os
import pandas as pd
from pyfasta import Fasta


class OrthogroupToGeneName():
    @staticmethod
    def run(path_to_orthogroups_tsv, path_to_fasta_dir, write=False):
        assert os.path.isfile(path_to_orthogroups_tsv), 'file --path_to_orthogroups_tsv does not exist: {}'.format(
            path_to_orthogroups_tsv)
        assert os.path.isdir(path_to_fasta_dir), 'folder --path_to_fasta_dir does not exist: {}'.format(
            path_to_fasta_dir)
        if write:
            out = os.path.dirname(path_to_orthogroups_tsv)
            assert os.path.isdir(out)

        # read Orthogroups.tsv
        df = pd.read_csv(path_to_orthogroups_tsv, sep='\t')
        df.set_index('Orthogroup', inplace=True)

        def get_gene_name(identifer: str):
            identifer = identifer.split(' ', maxsplit=1)[1].split(' [', maxsplit=1)[0]
            if identifer.startswith('hypothetical protein'):
                return 'hypothetical protein'
            else:
                return identifer

        def get_gene_id(identifer: str):
            return identifer.split(' ', maxsplit=1)[0]

        def get_gene_id_to_name_dict(strain):
            fasta_file_path = path_to_fasta_dir + '/{}.faa'.format(strain)
            assert os.path.isfile(fasta_file_path), 'fasta file {} is missing!'.format(fasta_file_path)
            gene_ids = [key for key in Fasta(fasta_file_path).keys()]
            return {get_gene_id(gene_id): get_gene_name(gene_id) for gene_id in gene_ids}

        def gene_id_to_gene_name(gene_ids, gene_id_to_name):
            if isinstance(gene_ids, float): return []
            return [gene_id_to_name[gene_id] for gene_id in gene_ids.split(', ')]

        strains = df.columns
        for strain in strains:
            gene_id_to_name = get_gene_id_to_name_dict(strain)
            df[strain] = df[strain].apply(gene_id_to_gene_name, args=([gene_id_to_name]))

        def majority_vote(row, best_only=False):
            all_names = [name for cell in row for name in cell]
            if best_only:
                return max(set(all_names), key=all_names.count)

            names_set = set(all_names)
            ddd = {gene_id: all_names.count(gene_id) for gene_id in names_set}
            return ddd.__str__()

        df_majority = pd.DataFrame(index=df.index)
        df_majority['Best Gene Name'] = df.apply(majority_vote, axis=1, args=[True])
        df_majority['Gene Name Occurrences'] = df.apply(majority_vote, axis=1)

        if write:
            out_path = out + '/Orthogroup_Best_Names.tsv'
            df_majority.to_csv(path_or_buf=out_path, sep='\t')
            print('Successfully wrote "{}"'.format(out_path))

        majority_dict = {orthogroup: majority_name for orthogroup, majority_name in
                         df_majority['Best Gene Name'].iteritems()}
        return majority_dict


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=
        """
        Create new file in the same folder as Orthogroups.tsv which contains 
        the most common gene name for each orthogroup.
        
        Usage as class:
        >>> from orthogroup_to_gene_name import OrthogroupToGeneName
        >>> OrthogroupToGeneName().run(
                path_to_orthogroups_tsv=<og_tsv>,
                path_to_fasta_dir=<fasta_dir>
                write=<True/False>
                )
            # returns dictionary: {'OG0000000': 'best-gene-name', ...}
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--og_tsv", type=str,
        help="Path to Orthogroups.tsv, usually located here: OrthoFinder/Results_{date}/Orthogroups/Orthogroups.tsv",
        required=True
    )
    parser.add_argument(
        "--fasta_dir", type=str, help="Path to folder where FASTA-files are stored. Fasta-files must end with .faa",
        required=True
    )
    args = parser.parse_args()

    OrthogroupToGeneName().run(path_to_orthogroups_tsv=args.og_tsv, path_to_fasta_dir=args.fasta_dir, write=True)
