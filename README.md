# get_tax_info

## Idea

* Calculate the most common gene name of an orthogroup by majority vote: `orthogroup_to_gene_name.py`
* Create plots analogous to roary_plots: `orthofinder_plots.py`

## Setup
Install dependencies, using your linux package manager or pip:
* orthogroup_to_gene_name.py: `pandas pyfasta fire`
* orthofinder_plots.py: `pandas argparse numpy biopython matplotlib seaborn`

## Usage
### orthogroup_to_gene_name.py
```
# Command line usage:

python3 orthogroup_to_gene_name.py --og_tsv path/to/Orthogroups.tsv --fasta_dir /path/to/fastas --write True
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
majority_dict = OrthogroupToGeneName.run(orthogroups_tsv=orthogroups_tsv, fasta_dir=fasta_dir, write=False)
```
`majority_dict` will be a python dict with key='orthogroup' -> value='best name', for example:

```
{
    'OG0000000': 'amino acid ABC transporter',
    'OG0000001': 'IS30 family transposase',
    'OG0000002': 'IS5/IS1182 family transposase',
}
```

If `write=True`, the script will create the `Orthogroup_BestNames.tsv` file as above and nothing will be returned.


### orthofinder_plots.py
**Disclaimer:**
This script is a port of [roary_plots](https://github.com/sanger-pathogens/Roary/tree/master/contrib/roary_plots) by Marco Galardini (marco@ebi.ac.uk).

```
# Command line usage:
python3 orthofinder_plots.py --tree data/SpeciesTree_rooted.txt --orthogroups_tsv data/Orthogroups.tsv --out=output
```

Three files will be created:

<img src="output/pangenome_frequency.svg"  width="80%"><br>
<img src="output/pangenome_matrix.svg"  width="80%"><br>
<img src="output/pangenome_pie.svg"  width="80%"><br>


#### Usage as python class
```
# load class
from orthofinder_plots import OrthofinderPlots

OrthofinderPlots.create_plots(
    phylo_object='/path/to/SpeciesTree_rooted.txt',
    pandas_table='/path/to/Orthogroups.tsv',
    output_format='svg',
    no_labels='False',
    output_location='/path/to/output/folder'
)
```
