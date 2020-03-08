import pickle
import re
import nltk 
from gensim.models import word2vec, Word2Vec
from gensim.models.phrases import Phrases, Phraser

def is_number(s):
    '''tests if a string can be converted in a number
    '''
    try:
        float(s)
        return True

    except ValueError:
        return False


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
    for k,ab in enumerate(pmid_Abstracts):
        key = ab['pmid']
        text = list(ab['title']) + [' '] + list(ab['abstract'])
        #check if abstract is empty 
        if len(ab['abstract']) == 0:
            continue
        
        #check if there are any known entities
        if 'entities' not in ab.keys():
            continue
            
        entities = ab['entities']
        for entity in entities:
            #check if entity list is empty
            if len(entities[entity])==0:
                continue
            
            ners = entities[entity]
            for ner in ners:
                #some objects are concept unique identifiers-less
                if ner['id']=='CUI-less':
                    continue
                
                for i in range(ner['start'], ner['end']):
                    text[i] = ' '
                    # +1 from start so always a whitespace before
                text[(ner['start']+1)] = str(ner['id'])
                
        #remove empty whitespace from filling in NER-IDs
        corpus[key] = re.sub('\s+', ' ', ''.join(text)).strip()
            
    return corpus          


def preprocess(corpora):
    '''Takes BERN token sutstituted text and does a final round 
        preprocessing prior to Word2Vec. Removes generic stop words 
        from NTLK (except for b and t, heurstic chosen). 

    Input:
    -corpora - dict, DOI's are key, entry is a string
    
    Return:
    -corpora_final - dict, but preprocessed with the above instructions
    '''
    corpora_final = {}
    stop_words = nltk.corpus.stopwords.words('english_biomed')
    
    for k, doi in enumerate(corpora):
        corpora[doi] = corpora[doi].lower()
        #remove everything except numbers and standard latin alphabet
        corpora[doi] = re.sub(r'[^a-zA-Z0-9\s]', '', corpora[doi])
        tokenize = corpora[doi].split(' ')
        tokenize = list(filter(None, tokenize))

        #UPDATE THIS USING str.endswith() method
        for ind, token in enumerate(tokenize):
            if token[-1]=='-':
                tokenize[ind] = token.split('-')[0]
            if token[0]=='-':
                tokenize[ind] = token.split('-')[1]
                
        #remove numerical values less than 8000 (lowest bernID = 8001)
            if is_number(token) is True:
                if float(token) < 8001: 
                    tokenize.pop(ind)

        #filter tokens to remove stopwords
        filtered_tokens = [token for token in tokenize if not token in stop_words]
        filtered_tokens = list(filter(None,filtered_tokens))
        
        # remove small corpuses from dict
        if len(filtered_tokens) > 24: 
            corpora_final[doi] = filtered_tokens
        
    return corpora_final

def phrase_corpus(corpus, min_count=5):
    '''Constructs phrases based from a certain threshold.
    Builds trigrams.
    
    Input: 
    - corpus: a list, each entry is an abstract

    Returns:
    - corpus: a list, with phrases in entries if applicable
    '''

    bigram = Phrases(corpus, min_count=min_count)
    trigram = Phrases(bigram[corpus], min_count=min_count)
    corpus = [trigram[bigram[corpus[k]]] for k, sent in enumerate(corpus)]
    return corpus


def train_word2vec(corpus, feature_size=200, 
                    window_context = 5, min_word_count=5,
                    iterations=1, sg=1):
    '''run word2vec on corpus, skip-gram setting automated
    default settings'''

    w2v_model = word2vec.Word2Vec(corpus, size=feature_size, 
                        window=window_context, min_count=min_word_count,
                        iter=iterations, sg=sg)

    return w2v_model

def main():

    # load data
    file = open('getPMID_Abstracts_out.txt', 'rb')
    corpus = pickle.load(file)

    corpus = process_PMIDabstracts(corpus)
    corpus = preprocess(corpus)

    # corpus is a dict, convert to list
    corpus = [corpus[key] for key in corpus]

    # build phrases
    corpus = phrase_corpus(corpus)†

    # generate word2vec and save
    w2v = train_word2vec(corpus)
    w2v.save("w2v_pbmc.model")

if __name__ == '__main__':
    main()