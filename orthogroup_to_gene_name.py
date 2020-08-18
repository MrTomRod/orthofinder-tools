import os
import logging
import pandas as pd
from Bio import SeqIO
from collections import Counter


class OrthogroupToGeneName:
    def __init__(self, n0_tsv: str, index_column='HOG'):
        assert os.path.isfile(n0_tsv), F'orthogroups_tsv does not exist: "{n0_tsv}"'
        # read "N0.tsv" and drop extra columns
        self.orig_df = pd.read_csv(n0_tsv, sep='\t')

        if index_column == 'HOG':
            self.orig_df.set_index('HOG', inplace=True)
            self.orig_df.drop(columns=['OG', 'Gene Tree Parent Clade'], inplace=True)
        elif index_column == 'OG':
            self.orig_df.set_index('OG', inplace=True)
            self.orig_df.drop(columns=['HOG', 'Gene Tree Parent Clade'], inplace=True)
        else:
            raise AssertionError(F'index_column must either be HOG or OG, but is {index_column}')

        self.strains = self.orig_df.columns

        self.translated_df = None  # use run to generate
        self.majority_df = None  # use run to generate
        self.majority_dict = None  # use run to generate

    def run(self, fasta_dir: str, file_endings='fasta', out=None):
        assert os.path.isdir(fasta_dir), F'fasta_dir does not exist: "{fasta_dir}"'
        if out is not None:
            assert os.path.isdir(os.path.dirname(out))

        self.translated_df = self.orig_df.__deepcopy__()

        def get_gene_id_to_name_dict(strain):
            fasta_file_path = os.path.join(fasta_dir + F'/{strain}.{file_endings}')
            assert os.path.isfile(fasta_file_path), F'{file_endings} file "{fasta_file_path}" is missing!'
            genes = SeqIO.parse(fasta_file_path, "fasta")
            return {gene.id: gene.description.split(' ', maxsplit=1)[1] for gene in genes}

        def gene_id_to_gene_name(gene_ids, gene_id_to_name):
            if isinstance(gene_ids, float): return []
            if gene_ids is None:
                return []
            try:
                return [gene_id_to_name[gene_id] for gene_id in gene_ids.split(', ')]
            except KeyError:
                logging.warning("Genes not present in fasta file: \n", gene_ids)
                return []

        for strain in self.strains:
            gene_id_to_name = get_gene_id_to_name_dict(strain)
            self.translated_df[strain] = self.translated_df[strain].apply(gene_id_to_gene_name, args=([gene_id_to_name]))

        def majority_vote(row):
            """ Get gene names set per HOG and count them. Disregards names with eAED in description, 
            typically useless maker annotations"""

            all_names = set(name for cell in row for name in cell if "eAED" not in name)
            if not all_names:
                all_names = ["NO DESCRIPTION"]  # if all names contain "eAED", all_names is an empty list

            names_to_count = Counter(all_names)

            return names_to_count.most_common(1)[0][0], str(names_to_count)  # return best gene name and gene name occurrences as string

        self.majority_df = pd.DataFrame(index=self.translated_df.index)
        self.majority_df[['Best Gene Name', 'Gene Name Occurrences']] = self.translated_df.apply(lambda row: pd.Series(majority_vote(row)), axis=1)

        if out is not None:
            self.majority_df.to_csv(path_or_buf=out, sep='\t')
            print(F'Successfully wrote "{out}"')
            return

        self.majority_dict = {
            orthogroup: majority_name
            for orthogroup, majority_name
            in self.majority_df['Best Gene Name'].iteritems()
        }

        return self.majority_dict


if __name__ == "__main__":
    import fire  # pip install fire # automated argparse


    def bootstrap(n0_tsv: str, fasta_dir: str, out_path: str, index_column='HOG', file_endings='fasta'):
        majority_dict = OrthogroupToGeneName(n0_tsv, index_column).run(fasta_dir, file_endings, out_path)


    fire.Fire(bootstrap)
