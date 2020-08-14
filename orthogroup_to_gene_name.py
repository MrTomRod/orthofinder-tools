import os
import pandas as pd
from Bio import SeqIO


class OrthogroupToGeneName():
    @staticmethod
    def run(n0_tsv: str, fasta_dir: str, file_endings='fasta', out_path=None, index_column='HOG'):
        assert os.path.isfile(n0_tsv), F'orthogroups_tsv does not exist: "{n0_tsv}"'
        assert os.path.isdir(fasta_dir), F'fasta_dir does not exist: "{fasta_dir}"'
        assert index_column in ['OG', 'HOG']

        if out_path is not None:
            assert os.path.isdir(os.path.dirname(out_path))

        def get_gene_name(identifer: str):
            identifer = identifer.split(' ', maxsplit=1)[1].split(' [', maxsplit=1)[0]
            if identifer.startswith('hypothetical protein'):
                return 'hypothetical protein'
            else:
                return identifer

        def get_gene_id_to_name_dict(strain):
            fasta_file_path = os.path.join(fasta_dir + F'/{strain}.{file_endings}')
            assert os.path.isfile(fasta_file_path), F'{file_endings} file "{fasta_file_path}" is missing!'
            genes = SeqIO.parse(fasta_file_path, "fasta")
            return {gene.id: gene.description.split(' ', maxsplit=1)[1] for gene in genes}

        def gene_id_to_gene_name(gene_ids, gene_id_to_name):
            if isinstance(gene_ids, float): return []
            if gene_ids is None: return []

            try:
                return [gene_id_to_name[gene_id] for gene_id in gene_ids.split(', ')]

            except KeyError:
                print("Genes not present in fasta file: \n", gene_ids)
                return []

        # read "N0.tsv" and drop extra columns
        df = pd.read_csv(n0_tsv, sep='\t')

        if index_column == 'HOG':
            df.set_index('HOG', inplace=True)
            df.drop(columns=['OG', 'Gene Tree Parent Clade'], inplace=True)
        elif index_column == 'OG':
            df.set_index('OG', inplace=True)
            df.drop(columns=['HOG', 'Gene Tree Parent Clade'], inplace=True)
        else:
            raise AssertionError(F'index_column must either be HOG or OG, but is {index_column}')

        strains = df.columns

        for strain in strains:
            gene_id_to_name = get_gene_id_to_name_dict(strain)
            df[strain] = df[strain].apply(gene_id_to_gene_name, args=([gene_id_to_name]))

        def majority_vote(row, best_only=False):
            """ Get gene names set per HOG and count them. Disregards names with eAED in description, 
            typically unusefull maker annotations"""

            all_names = [name for cell in row for name in cell if "eAED" not in name]
            if not all_names:  all_names = ["NO DESCRIPTION"]  # if all names contain "eAED", all_names is an empty list

            if best_only:
                return max(set(all_names), key=all_names.count)

            names_set = set(all_names)
            ddd = {gene_id: all_names.count(gene_id) for gene_id in names_set}
            return ddd.__str__()

        df_majority = pd.DataFrame(index=df.index)
        df_majority['Best Gene Name'] = df.apply(majority_vote, axis=1, args=[True])
        df_majority['Gene Name Occurrences'] = df.apply(majority_vote, axis=1)

        if out_path is not None:
            df_majority.to_csv(path_or_buf=out_path, sep='\t')
            print(F'Successfully wrote "{out_path}"')
            return

        majority_dict = {orthogroup: majority_name for orthogroup, majority_name in
                         df_majority['Best Gene Name'].iteritems()}
        return majority_dict


if __name__ == "__main__":
    import fire  # pip install fire # automated argparse

    fire.Fire(OrthogroupToGeneName.run)
