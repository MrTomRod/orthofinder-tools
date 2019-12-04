# get_tax_info

## Idea

* Calculate the most common gene name of an orthogroup by majority vote: `orthogroup_to_gene_name.py`
* Create plots analogous to roary_plots: `orthofinder_plots.py`

## Setup
Install dependencies, using your linux package manager or pip:
* orthogroup_to_gene_name.py: `pandas argparse pyfasta`
* orthofinder_plots.py: `pandas argparse numpy biopython matplotlib seaborn`

## Usage
### orthogroup_to_gene_name.py
```
# Command line usage:
python3 orthogroup_to_gene_name.py --og_tsv path/to/Orthogroups.tsv --fasta_dir /path/to/fastas
```

The product will be a file with the name `Orthogroup_BestNames.tsv` in the same folder as `Orthogroups.tsv`.

The tsv looks like this:

|   Orthogroup  |         Best Gene Name        | Gene Name Occurrences |
| ------------- | ----------------------------- | --------------------- |
| OG0000000     | amino acid ABC transporter    | {JSON}                |
| OG0000001     | IS30 family transposase       | {JSON}                |
| OG0000002     | IS5/IS1182 family transposase | {JSON}                |

The JSON is a dictionary with key='gene name' -> value=occurrence, for example:

```
{
    'Integrase core domain protein': 47,
    'hypothetical protein': 15,
    'IS30 family transposase': 126
}
```

#### Usage as python class
```
# load class
from orthogroup_to_gene_name import OrthogroupToGeneName
ogn = OrthogroupToGeneName()
majority_dict = ogn.run(path_to_orthogroups_tsv=args.og_tsv, path_to_fasta_dir=args.fasta_dir, write=False)
```

If `write=True`, the script will create the `Orthogroup_BestNames.tsv` file as above.

`majority_dict` will be a python dict with key='orthogroup' -> value='best name', for example:

```
{
    'OG0000000': 'amino acid ABC transporter',
    'OG0000001': 'IS30 family transposase',
    'OG0000002': 'IS5/IS1182 family transposase',
}

```

### orthofinder_plots.py
**Disclaimer:**
This script is a port of [roary_plots](https://github.com/sanger-pathogens/Roary/tree/master/contrib/roary_plots) by Marco Galardini (marco@ebi.ac.uk).

```
# Command line usage:
python3 orthofinder_plots.py --tree data/SpeciesTree_rooted.txt --spreadsheet data/Orthogroups.tsv --output_folder=output
```

Three files will be created:

![pangenome_frequency.svg](output/pangenome_frequency.svg)
![pangenome_matrix.svg](output/pangenome_matrix.svg)
![pangenome_pie.svg](output/pangenome_pie.svg)


#### Usage as python class
```
# load class
from orthofinder_plots import OrthofinderPlots
op = OrthofinderPlots()

tree = op.import_tree('path_to_tree')  # a phylo-object
pandas_table = op.import_orthofinder_table(args.spreadsheet)  # a pandas table

op.create_plots(
    phylo_object=tree,
    pandas_table=pandas_table,
    output_format='svg',
    no_labels='False',
    output_location='/path/to/output/folder'
)
```
