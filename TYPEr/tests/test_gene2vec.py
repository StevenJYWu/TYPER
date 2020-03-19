import gene2vec
import json
import pandas as pd

# load a dictionary that'll be used for testing
pmids = []
# load in small database of abstracts
file = open("test_db/pubmed19n1198.json", 'r')
data = [json.loads(line) for line in file]
for i in data:
    pmids.append(i['pmid'])


def test_is_number():
    """Test the function is_number, checks if a str can be represented as a
    float"""
    assert gene2vec.is_number("4324"), "float checker incorrect"


def test_getPMIDs():
    """Test whether or not we're correctly retreiving the pubmed ID numbers
    for abstracts based off of a search term we're interested in"""
    query = 'lol'
    nresults = 386
    results = gene2vec.getPMIDs(query)
    assert len(results) == nresults, "PMIDs not returned properly"


def test_getAbstracts():
    """Test whether or we're properly retrieving abstracts from a PMID search
    """
    corpus = gene2vec.getAbstracts(pmids, "test_db")
    # check if the correct number of abstracts is properly retrieved
    assert len(corpus) == 4935, "abstracts not properly retreived"


def test_process_PMIDabstracts():
    """Test to see integration of NER (named entity recognition) of words in
    our corpus is successful
    """
    string_match = '323201602 -Related 106929301'
    test_corp = [data[0]]
    corp = gene2vec.process_PMIDabstracts(test_corp)
    ls_strings = corp['31103025'].split(' ')
    first_3terms = ' '.join(ls_strings[0:3])
    assert first_3terms == string_match, "Abstracts processed incorrectly"


def test_preprocess():
    """Test if we're preprocessing our corpus correctly. Short abstracts are
    removed. See if it's functioning properly
    """
    corpus = {0: 'For writers, a random sentence can help them get their \
                creative juices flowing. Since the topic of the sentence is \
                completely unknown, it forces the writer to be creative when \
                the sentence appears. There are a number of different ways a \
                writer can use the random sentence for creativity.',
              1: 'THIS is a short sentence.'}

    corpus_pre = gene2vec.preprocess(corpus)
    assert len(corpus_pre[0]) == 26, \
        "Corpus tokenized/preprocessed incorrectly"
    try:
        corpus_pre[1]
        assert False, "short abstracts less than 25 words not removed properly"
    except KeyError:
        assert True


def test_phrase_corpus():
    """Test if we're properly generating phrases
    """
    documents = ["the mayor of new york was there",
                 "machine learning can be useful sometimes",
                 "new york is great",
                 "new york is worst than california",
                 "new york"]
    sentence_stream = [doc.split(" ") for doc in documents]
    sentence_stream = gene2vec.phrase_corpus(sentence_stream,
                                             min_count=2, threshold=1)
    assert sentence_stream[-1][0] == "new_york", \
        "phraser isn't working properly"


def test_train_word2vec():
    """unit testing a package library makes 0 sense
    """
    pass


def test_bernID_to_genes():
    """Testing bernID_to_gene method. The function is written to handle
    one specific file so unit testing is kind of excessive. Nevertheless
    we can try
    """
    path = "test_db/gene_meta_190805.tsv"
    df = pd.read_csv(path, sep='\t', header=None, dtype=str)
    bernID = gene2vec.bernID_to_genes(path)
    try:
        bernID[str(df.iloc[0][0])]
    except KeyError:
        assert False, "bernIDs not loaded properly"


def test_drop_embeddings():
    """This function is dependent on a w2v model. It just queries the model
    and drops keys that aren't relevant in our case all terms that aren't
    genes. Thus, there is no reason to unit test this function
    """
    pass


def test_get_sim_scores():
    """This function is dependent on a w2v model. It just queries top_n word
    embeddings that are most similar our entry. In our case the use would be
    to search a key (cell_type) and it will return the cosine distance of that
    word embedding to all genes
    """
    pass
