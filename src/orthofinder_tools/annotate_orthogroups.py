import json
import os
import pandas as pd
from Bio import SeqIO
from collections import Counter

from .utils import load_og, load_hog


class OrthogroupToGeneName:
    def __init__(self, fasta_dir: str, file_endings='fasta'):
        fasta_dir = os.path.abspath(os.path.expanduser(fasta_dir))
        assert os.path.isdir(fasta_dir), F'fasta_dir does not exist: "{fasta_dir}"'
        self.fasta_dir = fasta_dir
        self.file_endings = file_endings
        self.majority_df = None

    @property
    def majority_dict(self):
        assert hasattr(self, 'majority_df'), F'Load Orthogroups.tsv or N0.tsv first!'
        return {
            orthogroup: majority_name
            for orthogroup, majority_name
            in self.majority_df['Best Gene Name'].items()
        }

    def save_majority_df(self, out_file: str, header: bool = True):
        """
        Writes the following file:

        HOG Best Gene Name  Gene Name Occurrences
        N0.HOG0000000   amino acid ABC transporter Counter({'amino acid ABC transporter': 43})
        ...
        """
        out_file = os.path.abspath(out_file)
        assert hasattr(self, 'majority_df'), F'Load Orthogroups.tsv or N0.tsv first!'
        self.majority_df.to_csv(path_or_buf=out_file, sep='\t', header=header)

    def save_orthogroup_to_gene_ids(self, out_file: str):
        """
        Writes the following file (no header):

        N0.HOG0000000   gene_1, gene_2
        N0.HOG0000001   gene_3, gene_4, gene_5
        ...
        """
        out_file = os.path.abspath(out_file)
        with open(out_file, 'w') as f:
            for orthogroup, row in self.gene_ids_df.iterrows():
                gene_ids = [gene_id for gene_ids in row for gene_id in gene_ids]
                gene_ids = ', '.join(gene_ids)
                f.write(F'{orthogroup}\t{gene_ids}\n')

    def save_orthogroup_to_best_name(self, out_file: str, header: bool = False):
        """
        Writes the following file:

        Gene            Annotation
        N0.HOG0000000	amino acid ABC transporter ATP-binding protein
        N0.HOG0000001	ATP-binding cassette domain-containing protein
        ...
        """
        out_file = os.path.abspath(out_file)
        self.majority_df.drop(
            columns=['Gene Name Occurrences'],
            inplace=False
        ).to_csv(
            path_or_buf=out_file,
            sep='\t',
            index_label=['Gene'],
            header=['Annotation']
        )

    def load_og(self, og_tsv: str):
        self.gene_ids_df = load_og(og_tsv, result_type='gene-list')
        self.__load_gene_names()

    def load_hog(self, hog_tsv: str):
        self.gene_ids_df = load_hog(hog_tsv, result_type='gene-list')
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

    @staticmethod
    def __majority_vote(row):
        """ Get gene names set per HOG and count them. Disregards names with eAED in description,
        typically useless maker annotations"""

        all_names = [name for cell in row for name in cell if "eAED" not in name]
        if len(all_names) == 0:
            all_names = ["NO DESCRIPTION"]  # if all names contain "eAED", all_names is an empty list

        names_to_count = Counter(all_names)

        return names_to_count.most_common(1)[0][0], json.dumps(names_to_count)  # return best gene name and gene name occurrences as string

    def __get_gene_id_to_name_dict(self, strain):
        fasta_file_path = os.path.join(self.fasta_dir + F'/{strain}.{self.file_endings}')
        assert os.path.isfile(fasta_file_path), F'{self.file_endings} file "{fasta_file_path}" is missing!'

        def extract_description(gene):
            description = gene.description.split(' ', maxsplit=1)
            assert len(description) == 2, \
                F'Failed to extract description for strain={strain}:\n' \
                F'gene.id={gene.id}\n' \
                F'gene.description={gene.description}\n' \
                F'description={description}'
            description = description[1]
            if description.endswith(']'):
                # remove species description
                return description.rsplit(' [', maxsplit=1)[0]
            else:
                return description

        with open(fasta_file_path) as f:
            genes = SeqIO.parse(f, "fasta")
            res = {gene.id: extract_description(gene) for gene in genes}

        return res


def cli(
        orthogroups_tsv: str,
        fasta_dir: str,
        out: str,
        hog: bool = True,
        file_endings='fasta',
        simple: bool = True,
        header: bool = False,
):
    """
    Calculate the most common gene name of each orthogroup by majority vote

    :param orthogroups_tsv: path to Orthogroups.tsv or N0.tsv
    :param fasta_dir: path to where the protein FASTAS that OrthoFinder used as input are stored
    :param out: path to output file
    :param hog: if True: expect hierarchical orthogroup file (e.g.N0.tsv), if False: expect Orthogroups.tsv
    :param file_endings: file endings of the FASTA files (e.g. "faa", "fasta", "FASTA", etc)
    :param simple: if True: output maps orthogenes to best names; if False: orthogenes are mapped
    to best names and majority dict
    :param header: if True: add header in output file
    """
    otg = OrthogroupToGeneName(fasta_dir=fasta_dir, file_endings=file_endings)
    if hog:
        otg.load_hog(hog_tsv=orthogroups_tsv)
    else:
        otg.load_og(og_tsv=orthogroups_tsv)

    if simple:
        otg.save_orthogroup_to_best_name(out, header=header)
    else:
        otg.save_majority_df(out, header=header)


def main():
    import fire

    fire.Fire(cli)


if __name__ == '__main__':
    main()
