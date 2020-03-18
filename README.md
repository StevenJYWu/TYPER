# TYPEr

* Clone this Repo
```
cd
git clone https://github.com/StevenJYWu/TYPER.git
```


```
├── Documentation
├── LICENSE
├── README.md
├── TYPEr
│   ├── __init__.py
│   ├── core.py
│   ├── gene2vec.py
│   ├── tests
│   │   └── test_gene2vec.py
│   └── utils.py
├── docs
│   ├── CHEME\ 545_546\ Technical\ Review.pdf
│   ├── Proof_of_Concept.ipynb
│   ├── component_specifications.md
│   └── finalpresentation.pdf
├── examples
│   └── Demo.ipynb
└── requirements.txt

```



======
Usage:
`python gene2vec.py: -h`

```usage: gene2vec.py [-h] -q QUERY -d PATH -g ID -k CT [CT ...]

A program to mine open source research literature from pubmed. Generic use is
to query a cell type and the program returns celltype-gene associations based
off of word embeddings. All output will be written to a result file. The raw
corpus (not preprocessed), trained word2vec model, and celltype-gene
similarity scores will be stored there.

optional arguments:
  -h, --help            show this help message and exit
  -q QUERY, --Query QUERY
                        General Tissue/CellType of Interest
  -d PATH, --Database PATH
                        Path to location of database that stores named entity
                        recognized articles. Reference:
                        https://github.com/dmis-lab/bern
  -g ID, --GeneIDs ID   Path of location to BERN-entity IDs for genes
  -k CT [CT ...], --Celltype CT [CT ...]
                        A list of Keys to look up in the word embedding model.
                        Similarity scores are generated between these keys and
                        every gene.
                        
```
                        
BRIEF DESCRIPTION

## Environment

## Installation

## Usage 

DOCUMENTATION/INFORMATION ABOUT THE PACKAGE/MODULES/FUNCTIONS
