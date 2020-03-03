# COMPONENT SPECIFICATIONS
---

## 1. Data Extraction 

    1. User specifies cell types of interest. These are used as a query term in pubmed. 
    2. We retrieve the PMID's of all abstracts that are relevant to that specific query term. 
    3. We then retrieve the actual abstracts and title and build a corpus.
---
## 2. Preprocess Data
    
    1. Normalize all biological terms across the corpus. Do this by utilizing a bio-trained named entity recognition algorithm called bioBERT. 
    2. Further refine abstracts by removing non-english character words, punctuations, and stopwords. Only retain alpha-numerical characters. 
---
## 3. Train Word Embeddings

    1. Utilize word2vec in the genism package to train word vectors. 
    2. Build phrases of common terms
    3. We'll leverage supervised-learning by utilizing the skip-gram model

---
## 4. Process Word Embeddings

    1. Retrieve only embeddings that contain genes in their phrases
    2. Generate similarity scores between cell_types and genes
---
## 5. Data Integration/Visualization

    1. Integrate our learned word embeddings with scRNA-seq data. 
    2. Test if wordvec similarity scores are able to help improve classification of single cell RNA-sequences.