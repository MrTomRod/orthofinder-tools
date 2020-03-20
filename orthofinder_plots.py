# Copyright (C) <2019> University of Bern - Interfaculty Bioinformatics Unit

# This program is heavily based on Marco Galardinis (mgala@bu.edu) roary_plots.
# https://github.com/sanger-pathogens/Roary/tree/master/contrib/roary_plots

# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

__author__ = "Thomas Roder"  # thomas.roder@bioinformatics.unibe.ch
__version__ = '0.1.0'

import os
import pandas as pd
import numpy as np
from Bio import Phylo

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


class OrthofinderPlots():
    @staticmethod
    def import_tree(path_to_newick):
        assert os.path.isfile(path_to_newick)
        return Phylo.read(path_to_newick, 'newick')

    @staticmethod
    def import_orthofinder_table(path_to_table):
        assert os.path.isfile(path_to_table)
        pandas_table = pd.read_table(path_to_table, sep='\t', engine='c')
        pandas_table.set_index('Orthogroup', inplace=True)
        pandas_table = pandas_table.notna()  # make binary
        return pandas_table

    @staticmethod
    def import_roary_table(path_to_table, skipped_columns=14):
        assert os.path.isfile(path_to_table)
        pandas_table = pd.read_csv(path_to_table, sep=',', low_memory=False)

        # Set index (group name)
        pandas_table.set_index('Gene', inplace=True)
        # Drop the other info columns

        pandas_table.drop(list(pandas_table.columns[:skipped_columns - 1]), axis=1, inplace=True)

        # Transform it in a presence/absence matrix (1/0)
        pandas_table.replace('.{2,100}', 1, regex=True, inplace=True)
        pandas_table.replace(np.nan, 0, regex=True, inplace=True)

        return pandas_table.transpose()

    @staticmethod
    def create_plots(tree, orthogroups_tsv, out, format='svg', no_labels=False):
        assert os.path.isdir(out)
        assert format in ['png', 'tiff', 'pdf', 'svg']

        if type(tree) == str:
            tree = OrthofinderPlots.import_tree(tree)
        assert type(tree) == Phylo.Newick.Tree

        if type(orthogroups_tsv) == str:
            orthogroups_tsv = OrthofinderPlots.import_orthofinder_table(orthogroups_tsv)
        assert type(orthogroups_tsv) == pd.DataFrame

        matplotlib.use('Agg')
        sns.set_style('white')

        # Max distance to create better plots
        mdist = max([tree.distance(tree.root, x) for x in tree.get_terminals()])

        # Sort the matrix by the sum of strains presence
        idx = orthogroups_tsv.sum(axis=1).sort_values(ascending=False).index
        orthofinder_sorted = orthogroups_tsv.loc[idx]

        ## Plot pangenome frequency
        plt.figure(figsize=(7, 5))

        number_of_strains = orthogroups_tsv.shape[1]

        plt.hist(
            x=orthogroups_tsv.sum(axis=1),
            bins=number_of_strains,
            histtype="stepfilled",
            alpha=.7
        )

        plt.xlabel('No. of genomes')
        plt.ylabel('No. of genes')

        sns.despine(left=True,
                    bottom=True)

        save_path = os.path.join(out, 'pangenome_frequency.' + format)
        plt.savefig(save_path, dpi=300)
        plt.clf()

        # Sort the matrix according to tip labels in the tree
        tree_tip_labels = tree.get_terminals()
        orthofinder_sorted = orthofinder_sorted[[x.name for x in tree_tip_labels]]

        ## Plot presence/absence matrix against the tree
        with sns.axes_style('whitegrid'):
            fig = plt.figure(figsize=(17, 10))

            ax1 = plt.subplot2grid((1, 40), (0, 10), colspan=30)
            a = ax1.matshow(
                Z=orthofinder_sorted.T,
                cmap=plt.cm.Blues,
                vmin=0, vmax=1,
                aspect='auto',
                interpolation='none',
            )
            ax1.set_yticks([])
            ax1.set_xticks([])
            ax1.axis('off')

            ax = fig.add_subplot(1, 2, 1)
            # matplotlib v1/2 workaround
            try:
                ax = plt.subplot2grid((1, 40), (0, 0), colspan=10, facecolor='white')
            except AttributeError:
                ax = plt.subplot2grid((1, 40), (0, 0), colspan=10, axisbg='white')

            fig.subplots_adjust(wspace=0, hspace=0)

            ax1.set_title('Orthofinder matrix\n(%d gene clusters)' % orthogroups_tsv.shape[0])

            if no_labels:
                Phylo.draw(
                    tree=tree, axes=ax,
                    show_confidence=False,
                    label_func=lambda x: None,
                    xticks=([],), yticks=([],),
                    ylabel=('',), xlabel=('',),
                    xlim=(-mdist * 0.1, mdist + mdist * 0.1),
                    axis=('off',),
                    title=('Tree\n(%d strains)' % orthogroups_tsv.shape[1],),
                    do_show=False,
                )
            else:
                fsize = 12 - 0.1 * orthogroups_tsv.shape[1]
                if fsize < 7:
                    fsize = 7
                with plt.rc_context({'font.size': fsize}):
                    Phylo.draw(tree, axes=ax,
                               show_confidence=False,
                               label_func=lambda x: str(x)[:10],
                               xticks=([],), yticks=([],),
                               ylabel=('',), xlabel=('',),
                               xlim=(-mdist * 0.1, mdist + mdist * 0.45 - mdist * orthogroups_tsv.shape[1] * 0.001),
                               axis=('off',),
                               title=('Tree\n(%d strains)' % orthogroups_tsv.shape[1],),
                               do_show=False,
                               )

            save_path = os.path.join(out, 'pangenome_matrix.' + format)
            plt.savefig(save_path, dpi=300)
            plt.clf()

        ## Plot the pangenome pie chart
        plt.figure(figsize=(10, 10))

        core = orthogroups_tsv[(orthogroups_tsv.sum(axis=1) >= orthogroups_tsv.shape[1] * 0.99) & (
                orthogroups_tsv.sum(axis=1) <= orthogroups_tsv.shape[1])].shape[0]
        softcore = \
            orthogroups_tsv[(orthogroups_tsv.sum(axis=1) >= orthogroups_tsv.shape[1] * 0.95) & (
                    orthogroups_tsv.sum(axis=1) < orthogroups_tsv.shape[1] * 0.99)].shape[0]
        shell = orthogroups_tsv[(orthogroups_tsv.sum(axis=1) >= orthogroups_tsv.shape[1] * 0.15) & (
                orthogroups_tsv.sum(axis=1) < orthogroups_tsv.shape[1] * 0.95)].shape[
            0]
        cloud = orthogroups_tsv[orthogroups_tsv.sum(axis=1) < orthogroups_tsv.shape[1] * 0.15].shape[0]

        total = orthogroups_tsv.shape[0]

        def my_autopct(pct):
            val = int(round(pct * total / 100.0))
            return '{v:d}'.format(v=val)

        a = plt.pie(
            x=[core, softcore, shell, cloud],
            labels=['core\n(%d <= strains <= %d)' % (orthogroups_tsv.shape[1] * .99, orthogroups_tsv.shape[1]),
                    'soft-core\n(%d <= strains < %d)' % (
                        orthogroups_tsv.shape[1] * .95, orthogroups_tsv.shape[1] * .99),
                    'shell\n(%d <= strains < %d)' % (orthogroups_tsv.shape[1] * .15, orthogroups_tsv.shape[1] * .95),
                    'cloud\n(strains < %d)' % (orthogroups_tsv.shape[1] * .15)],
            explode=[0.1, 0.05, 0.02, 0], radius=0.9,
            colors=[(0, 0, 1, float(x) / total) for x in (core, softcore, shell, cloud)],
            autopct=my_autopct
        )

        save_path = os.path.join(out, 'pangenome_pie.' + format)
        plt.savefig(save_path, dpi=300)
        plt.clf()


if __name__ == "__main__":
    import fire
    fire.Fire(OrthofinderPlots.create_plots)
