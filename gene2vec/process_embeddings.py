from gensim.models import word2vec, Word2Vec
import pandas as pd
import numpy as np
import itertools
import pickle


def is_number(s):
    '''tests if a string can be converted in a number
    '''
    try:
        float(s)
        return True

    except ValueError:
        return False


def bernID_to_genes(df):
	'''This method is extremely hardcoded.
	Soley handles the file format of gene_extids_190510.tsv
	Inputs:
	- df: a dataframe, loaded from 
		'/Users/sjwu/nlp_CellMarkers/BERN_entity_IDs/gene_meta_190805.tsv'
	
	Return:
	- bernIDs: a dictionary, key is bern ID, entry is a tuple 
		(ensembl_id, gene)
	'''
	# string to search if ensembl ID exists for bernID
	ensembl = 'Ensembl'
	bernIDs = {}

	for k, i in df.iterrows():
		if ensembl in str(i[1]):
			key = i[3]
			ensembl_id = i[1].split(':')[-1]
			gene_id = i[2].split('|')[0]

			# remove gene_IDs with less than 2 characters
			if len(gene_id) < 3:
				continue

			bernIDs[key] = (ensembl_id, gene_id)

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
	lambda_filter = lambda x: (x[0] is not None) and\
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
		tup_tmp = bernID[key] 
		scores.append(tuple(tup_tmp[0], tup_tmp[1], score))

	return scores


def main():
	# load model trained from wordembeddings.py
	w2v_model = Word2Vec.load('w2v_pbmc.model')

	# load dataframe for bernID_to_genes method
	gene_key = pd.read_csv('BERN_entity_IDs/gene_extids_190510.tsv', 
							sep='\t', header=None)
	bernID = bernID_to_genes(gene_key)

	# key words to look up in w2v model 
	cell_types = ['cd34', 'cd14_monocytes', 'nk_cells', 't_reg', 'bcells']

	# update w2v_model
	w2v_model = drop_embeddings(w2v_model, cell_types, bernID)

	# get similarity scores
	itermodel = itertools.repeat(w2v_model, len(cell_types))
	iterct = itertools.repeat(cell_types, len(cell_types))
	iterIDs = itertools.repeat(bernID, len(cell_types))

	# map get_sim_scores function for each cell_type
	scores = list(map(get_sim_scores, itermodel, cell_types, iterct, iterIDs))
	scores = dict(zip(cell_types, scores))

	with open('scores.data', 'rb') as filehandle:
		pickle.dump(scores, filehandle)


	
if __name__ == '__main__':
	main()