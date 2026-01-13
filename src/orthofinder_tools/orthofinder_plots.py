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
__version__ = '0.2.0'

import json
import os, sys
import pandas as pd
import numpy as np
from Bio import Phylo

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from .utils import load_og, load_hog


def import_tree(path_to_newick):
    assert os.path.isfile(path_to_newick), f'File not found: {path_to_newick}'
    return Phylo.read(path_to_newick, 'newick')


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


def create_plots(tree, orthogroups_tsv, out, format='svg', no_labels=False, hog=False):
    """
    Create plots analogous to roary_plots for OrthoFinder OG/HOG

    :param tree: newick tree as string or path to file
    :param orthogroups_tsv: path to Orthogroups.tsv or N0.tsv
    :param out: path to output directory
    :param format: desired image format
    :param no_labels: Hide labels on phylogenetic tree
    :param hog: if True: expect hierarchical orthogroup file (e.g.N0.tsv), if False: expect Orthogroups.tsv
    """
    out = os.path.abspath(os.path.expanduser(out))
    os.makedirs(out, exist_ok=True)
    assert format in ['png', 'tiff', 'pdf', 'svg']

    if type(tree) == str:
        tree = os.path.abspath(tree)
        tree = import_tree(tree)
    assert type(tree) == Phylo.Newick.Tree

    if type(orthogroups_tsv) is str:
        if hog:
            orthogroups_tsv = load_hog(orthogroups_tsv, result_type='boolean')
        else:
            orthogroups_tsv = load_og(orthogroups_tsv, result_type='boolean')
    assert type(orthogroups_tsv) == pd.DataFrame

    matplotlib.use('Agg')
    sns.set_style('white')

    # Max distance to create better plots
    mdist = max([tree.distance(tree.root, x) for x in tree.get_terminals()])

    # Sort the matrix by the sum of strains presence
    og_count = orthogroups_tsv.sum(axis=1)
    idx = og_count.sort_values(ascending=False).index
    orthofinder_sorted = orthogroups_tsv.loc[idx]

    ## Plot pangenome frequency
    plt.figure(figsize=(7, 5))

    n_orthogenes, n_strains = orthogroups_tsv.shape

    plt.hist(
        x=og_count,
        bins=n_strains,
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

        ax = plt.subplot2grid((1, 40), (0, 0), colspan=10, facecolor='white')

        fig.subplots_adjust(wspace=0, hspace=0)

        ax1.set_title('OrthoFinder matrix\n(%d gene clusters)' % n_orthogenes)

        def draw_tree(xlim, label_func):
            Phylo.draw(
                tree=tree, axes=ax,
                show_confidence=False,
                label_func=label_func,
                xticks=([],), yticks=([],),
                ylabel=('',), xlabel=('',),
                xlim=xlim,
                axis=('off',),
                title=('Tree\n(%d strains)' % n_strains,),
                do_show=False,
            )

        if no_labels:
            draw_tree(
                xlim=(-mdist * 0.1, mdist + mdist * 0.1),
                label_func=lambda x: None
            )
        else:
            fsize = 12 - 0.1 * n_strains
            if fsize < 7:
                fsize = 7
            with plt.rc_context({'font.size': fsize}):
                draw_tree(
                    xlim=(-mdist * 0.1, mdist + mdist * 0.45 - mdist * n_strains * 0.001),
                    label_func=lambda x: str(x)[:10]
                )

        save_path = os.path.join(out, 'pangenome_matrix.' + format)
        plt.savefig(save_path, dpi=300)
        plt.clf()

    ## Plot the pangenome pie chart
    plt.figure(figsize=(10, 10))

    CORE, SOFT, SHELL = (n_strains * f for f in [.99, .95, .15])

    core = ((og_count >= CORE) & (og_count <= n_strains)).sum()
    softcore = ((og_count >= SOFT) & (og_count < CORE)).sum()
    shell = ((og_count >= SHELL) & (og_count < SOFT)).sum()
    cloud = (og_count < SHELL).sum()

    def my_autopct(pct):
        val = int(round(pct * n_orthogenes / 100.0))
        return '{v:d}'.format(v=val)

    pie_data = dict(zip(['core', 'softcore', 'shell', 'cloud'], [int(i) for i in (core, softcore, shell, cloud)]))
    print('Pie data:', json.dumps(pie_data))

    CORE, SOFT, SHELL = (int(i) for i in [CORE, SOFT, SHELL])
    ax = plt.pie(
        x=[core, softcore, shell, cloud],
        labels=[f'core\n({CORE} <= strains <= {n_strains})',
                f'soft-core\n({SOFT} <= strains < {CORE})',
                f'shell\n({SHELL} <= strains < {SOFT})',
                f'cloud\n(strains < {SHELL})'],
        explode=[.1, .05, .02, 0], radius=.9,
        colors=[(0, 0, 1, float(x) / n_orthogenes) for x in (core, softcore, shell, cloud)],
        autopct=my_autopct
    )

    save_path = os.path.join(out, 'pangenome_pie.' + format)
    plt.savefig(save_path, dpi=300)
    plt.clf()


def main():
    import fire

    fire.Fire(create_plots)


if __name__ == '__main__':
    main()
