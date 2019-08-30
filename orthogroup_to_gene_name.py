import os
import pandas as pd
from pyfasta import Fasta


class OrthogroupToGeneName():
    @staticmethod
    def run(path_to_orthogroups_tsv, path_to_fastas_dir, write=False):
        assert os.path.isfile(path_to_orthogroups_tsv)
        assert os.path.isdir(path_to_fastas_dir)
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
            gene_ids = [key for key in Fasta(path_to_fastas_dir + '/{}.faa'.format(strain)).keys()]
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
            df_majority.to_csv(path_or_buf=out + '/Orthogroup_Best_Names.tsv', sep='\t')

        majority_dict = {orthogroup: majority_name for orthogroup, majority_name in df_majority['Best Gene Name'].iteritems()}
        return majority_dict


path_to_orthogroups_tsv = '/home/thomas/PycharmProjects/44Paracasei/Genomics/faa/OrthoFinder/Results_Aug23_1/Orthogroups/Orthogroups.tsv'
path_to_fastas_dir = '/home/thomas/PycharmProjects/44Paracasei/Genomics/faa'

OrthogroupToGeneName().run(path_to_orthogroups_tsv, path_to_fastas_dir, write=True)
