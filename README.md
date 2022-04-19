# OrthofinderTools

## Idea

* Calculate the most common gene name of each orthogroup by majority vote: `orthogroup_to_gene_name.py`
* Create plots analogous to roary_plots: `orthofinder_plots.py`

## Setup
Install dependencies, using your linux package manager or pip:
* orthogroup_to_gene_name.py: `pandas biopython fire`
* orthofinder_plots.py: `pandas argparse numpy biopython matplotlib seaborn`

## Usage

### orthogroup_to_gene_name.py

#### Prerequisites

Your FASTA sequences must have some description, e.g.:

```text
>gnl|extdb|STRAIN-XY_000001 DNA-directed RNA polymerase subunit beta' [Pediococcus stilesii]
MIDVNKFESMQIGLASPDKIRMWSYGEVKKPETINYRTLKPEKDGLFDERIFGPTKDYECACGKYKRIRY
...
```

From this protein, `DNA-directed RNA polymerase subunit beta` will be extracted.

#### Command line usage

```
python3 orthogroup_to_gene_name.py \
    --orthofinder_tsv /path/to/N0_or_Orthogroups.tsv \
    --n0=True \  # True for N0.tsv, False for Orthogroups.tsv
    --fasta_dir /path/to/fastas \
    --out=outfile.tsv \
    --file_endings=faa \  # default=fasta; file suffix of the files in fasta_dir
```

The resulting tsv looks like this:

|   HOG         |         Best Gene Name        | Gene Name Occurrences |
| ------------- | ----------------------------- | --------------------- |
| N0.HOG0000000 | amino acid ABC transporter    | {JSON}                |
| N0.HOG0000001 | IS30 family transposase       | {JSON}                |
| N0.HOG0000002 | IS5/IS1182 family transposase | {JSON}                |

The JSON is a dictionary with key='gene name' -> value=occurrence, for example:

```json5
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

otg = OrthogroupToGeneName(
    fasta_dir=PATH_TO_ORTHOFINDER_FASTAS,
    file_endings='faa',
)
otg.load_hog(
    n0_tsv=F'{PATH_TO_ORTHOFINDER_FASTAS}/OrthoFinder/{CURRENT_FOLDER}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv'
)
```
`otg.majority_dict` will be a python dict with key='orthogroup' -> value='best name', for example:

```json5
{
    'N0.HOG0000000': 'amino acid ABC transporter',
    'N0.HOG0000001': 'IS30 family transposase',
    'N0.HOG0000002': 'IS5/IS1182 family transposase',
}
```

`otg.save_majority_df(outfile='path/to/outfile.tsv)` writes the following file:
```text
HOG Best Gene Name  Gene Name Occurrences
N0.HOG0000000   amino acid ABC transporter Counter({'amino acid ABC transporter': 43})
...
```

`otg.save_orthogroup_to_gene_ids(outfile='path/to/outfile.tsv)` writes the following file (no header):
```text
N0.HOG0000000   gene_1  gene_2
N0.HOG0000001   gene_3  gene_4  gene_5
...
```

`otg.save_orthogroup_to_gene_ids(outfile='path/to/outfile.tsv)` writes the following file (no header):
```text
N0.HOG0000000	amino acid ABC transporter ATP-binding protein
N0.HOG0000001	ATP-binding cassette domain-containing protein
...
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
```python
# load class
from orthofinder_plots import OrthofinderPlots

OrthofinderPlots.create_plots(
    tree='/path/to/SpeciesTree_rooted.txt',
    orthogroups_tsv='/path/to/Orthogroups.tsv',
    format='svg',
    no_labels=False,
    out='/path/to/output/folder'
)
```
