import os
import logging
import pandas as pd
from Bio import SeqIO
from collections import Counter


class OrthogroupToGeneName:
    def __init__(self, fasta_dir: str, file_endings='fasta'):
        assert os.path.isdir(fasta_dir), F'fasta_dir does not exist: "{fasta_dir}"'
        self.fasta_dir = fasta_dir
        self.file_endings = file_endings

    @property
    def majority_dict(self):
        assert hasattr(self, 'majority_df'), F'Load Orthogroups.tsv or N0.tsv first!'
        return {
            orthogroup: majority_name
            for orthogroup, majority_name
            in self.majority_df['Best Gene Name'].iteritems()
        }

    def save_majority_df(self, out_file: str):
        """
        Writes the following file:

        HOG Best Gene Name  Gene Name Occurrences
        N0.HOG0000000   amino acid ABC transporter Counter({'amino acid ABC transporter': 43})
        ...
        """
        assert hasattr(self, 'majority_df'), F'Load Orthogroups.tsv or N0.tsv first!'
        self.majority_df.to_csv(path_or_buf=out_file, sep='\t')

    def save_orthogroup_to_gene_ids(self, out_file: str):
        """
        Writes the following file (no header):

        N0.HOG0000000   gene_1, gene_2
        N0.HOG0000001   gene_3, gene_4, gene_5
        ...
        """
        with open(out_file, 'w') as f:
            for orthogroup, row in self.gene_ids_df.iterrows():
                gene_ids = [gene_id for gene_ids in row for gene_id in gene_ids]
                gene_ids = ', '.join(gene_ids)
                f.write(F'{orthogroup}\t{gene_ids}\n')

    def save_orthogroup_to_best_name(self, out_file: str):
        """
        Writes the following file (no header):

        N0.HOG0000000	amino acid ABC transporter ATP-binding protein
        N0.HOG0000001	ATP-binding cassette domain-containing protein
        ...
        """

        self.majority_df.drop(
            columns=['Gene Name Occurrences'],
            inplace=False
        ).to_csv(
            path_or_buf=out_file,
            sep='\t',
            header=False
        )

    def load_og(self, og_tsv: str):
        assert os.path.isfile(og_tsv), F'og_tsv does not exist: "{og_tsv}"'

        # read "Orthogroups.tsv"
        self.gene_ids_df = pd.read_csv(og_tsv, sep='\t', dtype=str)

        # sanity checks
        assert 'Orthogroup' in self.gene_ids_df.columns, \
            f'The file {og_tsv} does not appear to be a valid Orthogroup.tsv: ' \
            f'The column "Orthogroup" is missing. Try setting --n0=False'

        self.gene_ids_df.set_index('Orthogroup', inplace=True)
        self.gene_ids_df = self.gene_ids_df.applymap(lambda x: [] if pd.isnull(x) else x.split(', '))

        self.__load_gene_names()

    def load_hog(self, n0_tsv: str):
        assert os.path.isfile(n0_tsv), F'n0_tsv does not exist: "{n0_tsv}"'

        # read "N0.tsv" and drop extra columns
        self.gene_ids_df = pd.read_csv(n0_tsv, sep='\t', dtype=str)

        # sanity checks
        for col in ['HOG', 'OG', 'Gene Tree Parent Clade']:
            assert col in self.gene_ids_df.columns, f'The file {n0_tsv} does not appear to be a valid N0.tsv: ' \
                                                    f'The column "{col}" is missing. Try setting --n0=False'

        self.gene_ids_df.set_index('HOG', inplace=True)
        self.gene_ids_df.drop(columns=['OG', 'Gene Tree Parent Clade'], inplace=True)
        self.gene_ids_df = self.gene_ids_df.applymap(lambda x: [] if pd.isnull(x) else x.split(', '))

        self.__load_gene_names()

    def __load_gene_names(self):
        self.strains = self.gene_ids_df.columns

        self.gene_names_df = self.gene_ids_df.__deepcopy__()
        for strain in self.strains:
            gene_id_to_name = self.__get_gene_id_to_name_dict(strain)
            self.gene_names_df[strain] = self.gene_names_df[strain].apply(
                lambda ids: [gene_id_to_name[id] for id in ids])

        self.majority_df = pd.DataFrame(index=self.gene_names_df.index)
        self.majority_df[['Best Gene Name', 'Gene Name Occurrences']] = self.gene_names_df.apply(
            lambda row: pd.Series(self.__majority_vote(row)),
            axis=1)

    def __majority_vote(self, row):
        """ Get gene names set per HOG and count them. Disregards names with eAED in description,
        typically useless maker annotations"""

        all_names = [name for cell in row for name in cell if "eAED" not in name]
        if len(all_names) == 0:
            all_names = ["NO DESCRIPTION"]  # if all names contain "eAED", all_names is an empty list

        names_to_count = Counter(all_names)

        return names_to_count.most_common(1)[0][0], str(
            names_to_count)  # return best gene name and gene name occurrences as string

    def __get_gene_id_to_name_dict(self, strain):
        fasta_file_path = os.path.join(self.fasta_dir + F'/{strain}.{self.file_endings}')
        assert os.path.isfile(fasta_file_path), F'{self.file_endings} file "{fasta_file_path}" is missing!'
        genes = SeqIO.parse(fasta_file_path, "fasta")

        def extract_description(gene):
            description = gene.description.split(' ', maxsplit=1)
            assert len(description) == 2, \
                F'Failed to extract description for strain={strain}:\n' \
                F'gene.description={gene.description}\n' \
                F'description={description}'
            description = description[1]
            if description.endswith(']'):
                # remove species description
                return description.rsplit(' [', maxsplit=1)[0]
            else:
                return description

        return {gene.id: extract_description(gene) for gene in genes}


if __name__ == "__main__":
    import fire  # pip install fire # automated argparse


    def bootstrap(orthofinder_tsv: str, fasta_dir: str, out: str, n0: bool = True, file_endings='fasta'):
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

        otg.save_majority_df(out)


    fire.Fire(bootstrap)
