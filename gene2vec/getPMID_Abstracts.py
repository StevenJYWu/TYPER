
import requests
import re
from lxml import etree
import json
import glob
import pickle



def getPMIDs(query_term):
    """
    search a term in pubmed and get all PMIDs related to that query
    
    Input:
    -query_term: a string, to search PMIDs
    
    Output:
    -pmids: a list, of PMIDs
    """
    
    pmids = []
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term='
    #max number of PMIDs 
    params = {'retMax': 1000000} 
    #replace whitespaces with +
    query_term = re.sub('\s+', '+', query_term)
    
    response = requests.get(url=(url+query_term), params=params)
    if response.status_code == 200:
        root = etree.fromstring(response.content)
        IdList = root.find('IdList')
        for i in IdList:
            pmids.append(i.text)
    
    else:
        print('Requests has failed with error', response.status_code)
        
    return pmids


def getAbstracts(pmids):
    """
    Get BERN bioNER tagged abstracts from PMID
    Dependencies: The whole PUBMED database with BERN annotation downloaded from https://bern.korea.ac.kr/
        [to store in directory TYPEr/JSONS]
        
    *One caveat is that this database isn't up to date. Last updated in 08/2019
    
    Inputs:
    pmids: list, of pmids
    
    Outputs:
    corpus: list, of dictionaries
    """
    
    corpus = []
    jsons = glob.glob('[path to TYPEr/JSONS]')
    for k, js in enumerate(jsons):
        file = open(js, 'r')
        data = [json.loads(line) for line in file]
        for entry in data:
            if entry['pmid'] in pmids:
                corpus.append(entry)
                
        #progress check
        print(k)
    return corpus


search_term = 'peripheral blood mononuclear cells'
pmids = getPMIDs(search_term)
print("The numbe of Abstracts from search results are", len(pmids))
corpus = getAbstracts(pmids)

#save file as list
file = open("getPMID_Abstracts_out.txt", 'wb')   
pickle.dump(corpus, file)

