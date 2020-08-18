# OrthofinderTools

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

python3 orthogroup_to_gene_name.py \
    --n0_tsv /path/to/fastas/OrthoFinder/Results_Mon00/Phylogenetic_Hierarchical_Orthogroups/N0.tsv \
    --fasta_dir /path/to/fastas/OrthoFinder/fastas \
    --file_endings=fasta \   # sometimes faa
    --index_column=HOG \  # OG for orthogroups or HOG for more modern hierarchical orthogroups
    --out_path=/path/to/output.tsv
```

The tsv looks like this:

|   HOG         |         Best Gene Name        | Gene Name Occurrences |
| ------------- | ----------------------------- | --------------------- |
| N0.HOG0000000 | amino acid ABC transporter    | {JSON}                |
| N0.HOG0000001 | IS30 family transposase       | {JSON}                |
| N0.HOG0000002 | IS5/IS1182 family transposase | {JSON}                |

The JSON is a dictionary with key='gene name' -> value=occurrence, for example:

```
{
    'Integrase core domain protein': 47,
    'hypothetical protein': 15,
    'IS30 family transposase': 126
}
```

#### Usage as python class
```python
# load class
from orthogroup_to_gene_name import OrthogroupToGeneName

PATH_TO_ORTHOFINDER_FASTAS = '/path/to/OrthoFinder/fastas'
CURRENT_FOLDER = 'Results_Mon00'

majority_dict = OrthogroupToGeneName(
    n0_tsv=F'{PATH_TO_ORTHOFINDER_FASTAS}/OrthoFinder/{CURRENT_FOLDER}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv',
    index_column='OG'
).run(
    fasta_dir=PATH_TO_ORTHOFINDER_FASTAS,
    file_endings='faa',
)
```
`majority_dict` will be a python dict with key='orthogroup' -> value='best name', for example:

```
{
    'N0.HOG0000000': 'amino acid ABC transporter',
    'N0.HOG0000001': 'IS30 family transposase',
    'N0.HOG0000002': 'IS5/IS1182 family transposase',
}
```

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
