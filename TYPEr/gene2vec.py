import requests
import re
from lxml import etree
import json
import glob
import pickle
import nltk
from gensim.models import word2vec
from gensim.models.phrases import Phrases
import argparse
import os
import pandas as pd
import itertools

def is_number(s):
    '''tests if a string can be converted in a number
    '''
    try:
        float(s)
        return True

    except ValueError:
        return False


def getPMIDs(query_term):
    '''search a term in pubmed and get all PMIDs related to that query

    Input:
    -query_term: a string, to search PMIDs

    Output:
    -pmids: a list, of PMIDs
    '''

    pmids = []
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term='
    # max number of PMIDs
    params = {'retMax': 1000000}
    # replace whitespaces with +
    query_term = re.sub('\s+', '+', query_term)

    response = requests.get(url=(url + query_term), params=params)
    if response.status_code == 200:
        root = etree.fromstring(response.content)
        IdList = root.find('IdList')
        for i in IdList:
            pmids.append(i.text)

    else:
        print('Requests has failed with error', response.status_code)

    return pmids


def getAbstracts(pmids, dbpath):
    '''Get BERN bioNER tagged abstracts from PMID
    Dependencies: The whole PUBMED database with BERN annotation downloaded from https://bern.korea.ac.kr/
        -Stored on my local drive under: /Users/sjwu/nlp_CellMarkers/data/2019_merged_json_fixed

    *One caveat is that this database isn't up to date. Last updated in 08/2019

    Inputs:
    pmids: list, of pmids
    dbpath: string, path to annotated abstracts

    Return:
    corpus: list, of dictionaries'''
    # convert pmids to dict for faster lookup
    pmids = {i:'' for i in pmids}

    corpus = []
    jsons = glob.glob((dbpath + '*.json'))
    for k, js in enumerate(jsons):
        file = open(js, 'r')
        data = [json.loads(line) for line in file]

        for entry in data:
            if entry['pmid'] in pmids:
                corpus.append(entry)
    return corpus


def process_PMIDabstracts(pmid_Abstracts):
    '''converts PMID abstracts into BERN annotated abstracts.
        Replaces the strings with BioNER terms IDs.

        *Note: for getPMID_Abstracts.py, the dataebase reference
        sometimes have multiple occurances. Meaning PMID's/abstracts
        in the load_file are sometimes duplicated. Using a dictionary
        with PMID as a key, resolves this issue.

        *Note 2: The corpus object, includes title+abstracts in each entry

    Inputs:
    -pmid_Abstracts: list, of outputs from getPMID_Abstracts.py

    Return:
    -corpus: dictionary, PMID key as entry
    '''

    corpus = {}
    for k, ab in enumerate(pmid_Abstracts):
        key = ab['pmid']
        text = list(ab['title']) + [' '] + list(ab['abstract'])
        # check if abstract is empty
        if len(ab['abstract']) == 0:
            continue

        # check if there are any known entities
        if 'entities' not in ab.keys():
            continue

        entities = ab['entities']
        for entity in entities:
            # check if entity list is empty
            if len(entities[entity]) == 0:
                continue

            ners = entities[entity]
            for ner in ners:
                # some objects are concept unique identifiers-less
                if ner['id'] == 'CUI-less':
                    continue

                for i in range(ner['start'], ner['end']):
                    text[i] = ' '
                    # +1 from start so always a whitespace before
                text[(ner['start'] + 1)] = str(ner['id'])

        # remove empty whitespace from filling in NER-IDs
        corpus[key] = re.sub('\s+', ' ', ''.join(text)).strip()

    return corpus


def preprocess(corpora):
    '''Takes BERN token sutstituted text and does a final round
        preprocessing prior to Word2Vec. Removes generic stop words
        from NTLK.

    Input:
    -corpora - dict, DOI's are key, entry is a string

    Return:
    -corpora_final - dict, but preprocessed with the above instructions
    '''
    corpora_final = {}
    stop_words = nltk.corpus.stopwords.words('english')

    for k, doi in enumerate(corpora):
        corpora[doi] = corpora[doi].lower()
        # remove everything except numbers and standard latin alphabet
        corpora[doi] = re.sub(r'[^a-zA-Z0-9\s]', '', corpora[doi])
        tokenize = corpora[doi].split(' ')
        tokenize = list(filter(None, tokenize))
        # UPDATE THIS USING str.endswith() method
        for ind, token in enumerate(tokenize):
            if token[-1] == '-':
                tokenize[ind] = token.split('-')[0]
            if token[0] == '-':
                tokenize[ind] = token.split('-')[1]

            # remove numerical values less than 8000 (lowest bernID = 8001)
            if is_number(token) is True:
                if float(token) < 8001:
                    tokenize.pop(ind)

        # filter tokens to remove stopwords
        filtered_tokens = [token for token in tokenize if not token in stop_words]
        filtered_tokens = list(filter(None, filtered_tokens))
        # remove small corpuses from dict
        if len(tokenize) > 24:
            corpora_final[doi] = filtered_tokens

    return corpora_final


def phrase_corpus(corpus, min_count=5, threshold=5):
    '''Constructs phrases based from a certain threshold.
    Builds trigrams.

    Input:
    - corpus: a list, each entry is an abstract

    Returns:
    - corpus: a list, with phrases in entries if applicable
    '''

    bigram = Phrases(corpus, min_count=min_count, threshold=threshold)
    trigram = Phrases(bigram[corpus], min_count=min_count, threshold=threshold)
    corpus = [trigram[bigram[corpus[k]]] for k, sent in enumerate(corpus)]
    return corpus


def train_word2vec(corpus, feature_size=200,
                   window_context=5, min_word_count=5,
                   iterations=1, sg=1, negative=5):
    '''run word2vec on corpus, skip-gram setting automated
    default settings'''
    print(negative)
    w2v_model = word2vec.Word2Vec(corpus, size=feature_size,
                                  window=window_context, min_count=min_word_count,
                                  iter=iterations, sg=sg, negative=negative)
    return w2v_model


def bernID_to_genes(id_path):
    '''This method is extremely hardcoded.
    Soley handles the file format of gene_extids_190510.tsv (aka args.id)

    Inputs:
    - id_path: a string, path to gene_ID file

    Return:
    - bernIDs: a dictionary, key is bern ID, entry is a tuple
        (ensembl_id, gene)
    '''

    df = pd.read_csv(id_path, sep='\t', header=None, dtype=str)

    # string to search if ensembl ID exists for bernID
    ensembl = 'Ensembl'
    bernIDs = {}

    for k, i in df.iterrows():
        if ensembl in str(i[4]):
            key = i[0]
            ensembl_id = i[4].split(':')[-1]
            gene_id = i[3].split('|')[0]

            # remove gene_IDs with less than 2 characters
            if len(gene_id) < 3:
                continue

            bernIDs[str(key)] = (ensembl_id, gene_id)
        else:
            continue

    return bernIDs


def drop_embeddings(w2v_model, cell_types, bernID):
    ''' Drop all embeddings that don't have a celltype
    or gene
    Inputs:
    - w2v_model: a trained word2vec model
    - cell_types: a list of cell_type keys to retain
    - bernID: a dictionary, key is bern ID, entry is a tuple
        (ensembl_id, gene)

    Returns:
    - w2v_model: updated w2v_model with only useful terms '''

    # define all terms
    keys = list(w2v_model.wv.vocab.keys())

    # filter out non_gene terms
    get_genes = lambda x: any(str(i) in x for i in bernID)
    keys_gene = list(filter(get_genes, keys))

    # valid list of terms
    valid_keys = keys_gene + cell_types
    for idx, word in enumerate(w2v_model.wv.index2word):
        if word not in valid_keys:
            w2v_model.wv.index2word[idx] = None

    return w2v_model


def get_sim_scores(w2v_model, key, cell_types, bernID):
    ''' get similarity scores for a key filtering all
    none gene entries.

    Inputs:
    - w2v_model: a trained word2vec model
    - key: a string, of a vocab in the model
    - cell_types: a list of cell_type keys
    - bernID: a dictionary, key is bern ID, entry is a tuple
        (ensembl_id, gene)

    Returns:
    - scores: a list of (ensembl_id, gene_id, cos_score)
    '''
    nvocabs = len(w2v_model.wv.vocab)
    try:
        sim_scores = w2v_model.wv.most_similar(key, topn=nvocabs)

    # catch cases where key isn't in model
    except:
        print('key is not in model')
        return [None]

    # Filter out "None Values and other cell_types"
    lambda_filter = lambda x: (x[0] is not None) and \
                              (x[0] not in cell_types)
    scores = list(filter(lambda_filter, sim_scores))

    # break up phrases and assign the highest score to each phrase
    dict_tmp = {}
    for tup in scores:
        # tup has form ('phrase', cosine_score)
        words = tup[0].split('_')
        for word in words:
            isnumber = is_number(word)
            if (is_number(word) is True) and (word not in dict_tmp):
                dict_tmp[word] = tup[1]

    scores = []
    for key in dict_tmp:
        score = dict_tmp[key]
        if key in bernID:
            tup_tmp = bernID[key]
            scores.append((tup_tmp[0], tup_tmp[1], score))
        else:
            continue

    return scores


def main():
    pmids = getPMIDs(args.query)
    corpus = getAbstracts(pmids, args.path)

    # make results folder
    os.makedirs('result')

    # save file as list
    file = open("result/raw_corpus.bin", 'wb')
    pickle.dump(corpus, file)

    # preprocess corpus
    corpus = process_PMIDabstracts(corpus)
    corpus = preprocess(corpus)

    # corpus is a dict, convert to list
    corpus = [corpus[key] for key in corpus]

    # build phrases
    corpus = phrase_corpus(corpus, min_count=10, threshold=15)

    # generate word2vec and save
    w2v_model = train_word2vec(corpus, min_word_count=10, negative=15)

    bernID = bernID_to_genes(args.id)
    w2v_model = drop_embeddings(w2v_model, args.ct, bernID)
    w2v_model.save("result/w2v.model")

    # get similarity scores
    itermodel = itertools.repeat(w2v_model, len(args.ct))
    iterct = itertools.repeat(args.ct, len(args.ct))
    iterIDs = itertools.repeat(bernID, len(args.ct))

    # map get_sim_scores function for each cell_type
    scores = list(map(get_sim_scores, itermodel, args.ct, iterct, iterIDs))
    scores = dict(zip(args.ct, scores))

    with open('result/w2v_ct.scores', 'wb') as filehandle:
        pickle.dump(scores, filehandle)


parser = argparse.ArgumentParser(description='A program to mine open source\
            research literature from pubmed. Generic use is to query a cell\
            type and the program returns celltype-gene associations based\
            off of word embeddings. All output will be written to a result \
            file. The raw corpus (not preprocessed), trained word2vec model,\
            and celltype-gene similarity scores will be stored there.')
parser.add_argument('-q', '--Query', help='General Tissue/CellType of Interest',
                    dest='query', required=True)
parser.add_argument('-d', '--Database', help='Path to location of database that\
                    stores named entity recognized articles. Reference:\
                    https://github.com/dmis-lab/bern',
                    dest='path', required=True)
parser.add_argument('-g', '--GeneIDs', help='Path of location to BERN-entity \
                    IDs for genes', dest='id', required=True)
parser.add_argument('-k', '--Celltype', nargs='+',
        help='A list of Keys to look up in the word embedding model.\
             Similarity scores are generated between these keys and every\
             gene.', dest='ct', required=True)
args = parser.parse_args()

main()

