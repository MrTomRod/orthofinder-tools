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

import json
import os
import pandas as pd
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from .utils import load_og, load_hog


def import_tree(path_to_newick):
    if not os.path.isfile(path_to_newick):
        raise FileNotFoundError(f'Tree file not found: {path_to_newick}')
    try:
        return Phylo.read(path_to_newick, 'newick')
    except Exception as e:
        raise ValueError(f"Could not read tree file {path_to_newick}: {e}")


def import_roary_table(path_to_table, skipped_columns=14):
    if not os.path.isfile(path_to_table):
        raise FileNotFoundError(f"Roary table file not found: {path_to_table}")

    try:
        pandas_table = pd.read_csv(path_to_table, sep=',', low_memory=False)
    except Exception as e:
        raise ValueError(f"Could not read roary table {path_to_table}: {e}")

    if 'Gene' not in pandas_table.columns:
        raise ValueError(f"Roary table {path_to_table} is missing 'Gene' column")

    # Set index (group name)
    pandas_table.set_index('Gene', inplace=True)

    # Drop the other info columns
    cols_to_drop = list(pandas_table.columns[:skipped_columns - 1])
    pandas_table.drop(cols_to_drop, axis=1, inplace=True)

    # Transform it in a presence/absence matrix (1/0)
    pandas_table.replace('.{2,100}', 1, regex=True, inplace=True)
    pandas_table.replace(np.nan, 0, regex=True, inplace=True)

    return pandas_table.transpose()


def plot_pangenome_frequency(og_count, n_strains, out_path):
    """Plot pangenome frequency histogram."""
    fig = plt.figure(figsize=(7, 5))
    plt.hist(
        x=og_count,
        bins=n_strains,
        histtype="stepfilled",
        alpha=.7
    )
    plt.xlabel('No. of genomes')
    plt.ylabel('No. of genes')
    sns.despine(left=True, bottom=True)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_pangenome_matrix(tree, orthofinder_sorted, n_orthogenes, n_strains, out_path,
                          figsize=(17, 10), no_labels=False, tree_width=1, matrix_width=3, wspace=0.0):
    """Plot presence/absence matrix against the tree."""
    mdist = max([tree.distance(tree.root, x) for x in tree.get_terminals()])

    with sns.axes_style('whitegrid'):
        fig = plt.figure(figsize=figsize)

        # Use GridSpec with minimal margins
        gs = GridSpec(1, 2, width_ratios=[tree_width, matrix_width], wspace=wspace,
                      left=0.02, right=0.98, top=0.95, bottom=0.05)

        ax_matrix = fig.add_subplot(gs[1])
        ax_matrix.matshow(
            Z=orthofinder_sorted.T,
            cmap=plt.cm.Blues,
            vmin=0, vmax=1,
            aspect='auto',
            interpolation='none',
        )
        ax_matrix.set_yticks([])
        ax_matrix.set_xticks([])
        ax_matrix.axis('off')
        ax_matrix.set_title('OrthoFinder matrix\n(%d gene clusters)' % n_orthogenes)

        ax_tree = fig.add_subplot(gs[0], facecolor='white')

        def draw_tree(xlim, label_func):
            Phylo.draw(
                tree=tree, axes=ax_tree,
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
                xlim=(-mdist * 0.02, mdist + mdist * 0.02),
                label_func=lambda x: None
            )
        else:
            fsize = 12 - 0.1 * n_strains
            if fsize < 7:
                fsize = 7
            with plt.rc_context({'font.size': fsize}):
                draw_tree(
                    xlim=(-mdist * 0.02, mdist + mdist * 0.2),
                    label_func=str
                )

        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close(fig)


def plot_pangenome_pie(og_count, n_orthogenes, n_strains, out_path):
    """Plot the pangenome pie chart."""
    fig = plt.figure(figsize=(10, 10))

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

    CORE_int, SOFT_int, SHELL_int = (int(i) for i in [CORE, SOFT, SHELL])
    plt.pie(
        x=[core, softcore, shell, cloud],
        labels=[f'core\n({CORE_int} <= strains <= {n_strains})',
                f'soft-core\n({SOFT_int} <= strains < {CORE_int})',
                f'shell\n({SHELL_int} <= strains < {SOFT_int})',
                f'cloud\n(strains < {SHELL_int})'],
        explode=[.1, .05, .02, 0], radius=.9,
        colors=[(0, 0, 1, float(x) / n_orthogenes) for x in (core, softcore, shell, cloud)],
        autopct=my_autopct
    )

    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def create_plots(
        tree,
        orthogroups_tsv,
        out,
        format='svg',
        figsize="17,10",
        no_labels=False,
        hog=False,
        tree_width=1,
        matrix_width=3,
        wspace=0.0
):
    """
    Create plots analogous to roary_plots for OrthoFinder OG/HOG

    :param tree: newick tree as string or path to file
    :param orthogroups_tsv: path to Orthogroups.tsv or N0.tsv
    :param out: path to output directory
    :param format: desired image format
    :param figsize: size of the presence-absence-figure, e.g. "17,10"
    :param no_labels: Hide labels on phylogenetic tree
    :param hog: if True: expect hierarchical orthogroup file (e.g.N0.tsv), if False: expect Orthogroups.tsv
    :param tree_width: relative width of the tree plot
    :param matrix_width: relative width of the matrix plot
    :param wspace: space between the tree and the matrix
    """
    out = os.path.abspath(os.path.expanduser(out))
    os.makedirs(out, exist_ok=True)
    assert format in ['png', 'tiff', 'pdf', 'svg'], f'{format=} is invalid'

    if isinstance(tree, str):
        tree = import_tree(os.path.abspath(tree))

    if not isinstance(tree, Phylo.Newick.Tree):
        raise TypeError(f"tree must be a Newick Tree object or a path to a newick file, got {type(tree)}")

    if type(orthogroups_tsv) is str:
        if hog:
            orthogroups_tsv = load_hog(orthogroups_tsv, result_type='boolean')
        else:
            orthogroups_tsv = load_og(orthogroups_tsv, result_type='boolean')

    if not isinstance(orthogroups_tsv, pd.DataFrame):
        raise TypeError(
            f"orthogroups_tsv must be a pandas DataFrame or a path to a tsv file, got {type(orthogroups_tsv)}")

    if isinstance(figsize, str):
        try:
            figsize = tuple(float(x) for x in figsize.split(','))
        except ValueError:
            raise ValueError(f"figsize must be a comma-separated pair of numbers, got {figsize}")

    sns.set_style('white')

    # Data preparation
    n_orthogenes, n_strains = orthogroups_tsv.shape
    og_count = orthogroups_tsv.sum(axis=1)

    # Sort for matrix plot
    idx = og_count.sort_values(ascending=False).index
    orthofinder_sorted = orthogroups_tsv.loc[idx]
    tree_tip_labels = tree.get_terminals()
    orthofinder_sorted = orthofinder_sorted[[x.name for x in tree_tip_labels]]

    # 1. Pangenome frequency
    plot_pangenome_frequency(
        og_count, n_strains,
        os.path.join(out, 'pangenome_frequency.' + format)
    )

    # 2. Presence/absence matrix
    plot_pangenome_matrix(
        tree, orthofinder_sorted, n_orthogenes, n_strains,
        os.path.join(out, 'pangenome_matrix.' + format),
        figsize=figsize, no_labels=no_labels,
        tree_width=tree_width, matrix_width=matrix_width, wspace=wspace
    )

    # 3. Pangenome pie chart
    plot_pangenome_pie(
        og_count, n_orthogenes, n_strains,
        os.path.join(out, 'pangenome_pie.' + format)
    )


def main():
    import fire

    fire.Fire(create_plots)


if __name__ == '__main__':
    main()
