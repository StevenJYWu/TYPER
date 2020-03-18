![test image size](https://images-na.ssl-images-amazon.com/images/I/71JDzraloRL._AC_SY355_.jpg)

## TYPER 
This package's purpose is to generate celltype-gene association scores. The scores can then be uesd to do whatever downstream task the user is interested in. For example, under /examples we provide a demo.iypnb that takes you through one potential use of these score. In our demo, we utilize these associations to annotate single cell RNA-sequencing data. 


## Software Dependencies
* Python3
* See requirements.txt

## Data Dependencies 
* TYPER is dependent on an algorithm called BERN, it's a [Neural Named Entity Recognition and Multi-Type Normalization Tool for Biomedical Text Mining](https://bern.korea.ac.kr/). 
    
    * BERN is to computationally expensive to run on a local CPU, so for the purpose of this project we took advantage of the fact that they already applied their algorithm to all 20 million articles availabe in pubmed. The database reference is required to run our package and can be [downloaded here](https://drive.google.com/open?id=14YrlOGd1NdDn0XD-Yat4bbq3lRv1EyqR).
    
* BERN applies a NER algorithm and replaces redundent terms in the corpus with a concept unique identifier (CUI). A CUI-list annotation is necessary as well and can be [found here](https://drive.google.com/open?id=1KgJPBYB8D4_hN7wbiu0XOOM-lQdV8EgP).


## Organization of the Project
```
├── Documentation
├── LICENSE
├── README.md
├── TYPEr
│   ├── __init__.py
│   ├── core.py
│   ├── gene2vec.py
│   └── tests
│       └── test_gene2vec.py
├── docs
│   ├── CHEME\ 545_546\ Technical\ Review.pdf
│   ├── Proof_of_Concept.ipynb
│   ├── Untitled.ipynb
│   ├── component_specifications.md
│   └── finalpresentation.pdf
├── examples
│   ├── demo.ipynb
│   └── utils.py
└── requirements.txt

```

## Getting Started
* Clone this Repo
```
cd
git clone https://github.com/StevenJYWu/TYPER.git
```

* Install python packages
```
pip3 install -r requirements.txt --user
```

* Install Data Dependencies illustrated above

* Usage:

`python gene2vec.py: -h
`

```
usage: gene2vec.py [-h] -q QUERY -d PATH -g ID -k CT [CT ...]

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
