import pandas as pd
import numpy as np
from sklearn.decomposition import TruncatedSVD
import sklearn.feature_extraction.text as text
from sklearn.metrics import accuracy_score


def update_fdata(fdata, data):
    '''function to update fdata dataframe with NLP scores for each celltype

    Input:
    - fdata: a dataframe, with a column of valid ensembl ID's
    - data: dictionary[cell_type query] = list of (ensembl_id, gene_id, score)

    Returns:
    - fdata: updated dataframe, with additional columns of Celltype_nlp scores

    '''
    for key in data:
        list_tmp = []
        df = data[key]
        ref_ids = list(df.ensembl)

        for k, row in fdata.iterrows():
            ids = row[0]
            if ids in ref_ids:
                score = df[df['ensembl'] == ids].score
                list_tmp.append(float(score))
            else:
                list_tmp.append(0)

        fdata[key] = list_tmp

    return fdata


def format_data(data, fdata):
    '''function to format data loaded from gene2vec.py

    Input:
    - data: dictionary[cell_type query] = list of (ensembl_id, gene_id, score)
    - fdata: a dataframe, with a column of valid ensembl ID's

    Return:
    - data: updated dictionary with filtered (for valid ensembl ID)\
    and sorted dataframes
    '''

    for key in data:
        df = pd.DataFrame(data[key], columns=['ensembl', 'gene_id', 'score'])

        # filter out non_valid ensembl ID's
        df = df[df['ensembl'].isin(fdata['ensembl'])]
        df = df.sort_values('ensembl')
        data[key] = df
    return data


def downsample(pdata, N=2500):
    '''Down sample so that there is an equivalent number of cells from each type

    Inputs:
    - pdata: a dataframe, with a barcode and cell_type column
    - N: number of single cells from each type

    Returns:
    - idx: tuple of index of all selected subsamples
    '''
    idx = []
    for ct in pdata['celltype'].unique():
        ct_idx = pdata[pdata['celltype'] == ct].sample(N).index
        idx += list(ct_idx)

    return tuple(idx)


def l2_norm(X):
    '''Calculates the l2 norm across each row samples
    Inputs:
    - X: a matrix, shape (NxM)

    Return:
    - dist: a matrix, shape (NxN)
    '''
    N = X.shape[0]
    dist = np.zeros([N, N])
    for i in range(N):
        dist[i, :] = np.linalg.norm(X - X[i, :], axis=1)
    return dist


def most_frequent(List):
    '''Returns the most frequent occurence from a list'''
    return max(set(List), key=List.count)


def NN_search(graph, n, naive_ct):
    '''Search n-Nearest-Neighbor and vote on single cell identity
    Inputs:
    - graph: a NxN matrix depicting L2 \
    distance of each point from every other point
    - n: integer, number of neighbors to select
    - naive_ct: list, of naively defined celltypes

    Return:
    - nn_ct: a list of near-neighbor defined celltypes
    '''
    rows = graph.shape[0]
    nn_ct = []

    for i in range(rows):
        closestn = graph[i, :].argsort()[:n]
        closestn = list(naive_ct[closestn])
        nn_ct.append(most_frequent(closestn))
    return nn_ct


def tf_idf(matrix):
    '''converts a 2-d matrix to TF-IDF format

    Inputs:
    - matrix: 2-d array or sparse matrix

    Outputs:
    - tf_idf: 2-d array or sparse matrix
    '''
    transformer = text.TfidfTransformer(smooth_idf=True,
                                        use_idf=True, sublinear_tf=True)
    tfidf = transformer.fit_transform(matrix)
    return tfidf


def lsa(matrix, n_components=50):
    '''Returns a n_components-low rank matrix approximation, \
    purpose is to reduce noise

    Inputs:
    - matrix: 2-d array or sparse matrix

    Returns:
    - matrix: (N x n_components) shape matrix
    '''
    svd = TruncatedSVD(n_components=n_components, random_state=42)
    return svd.fit_transform(matrix.T)


def accuracy(y, scores, map_dict):
    '''Calculates accuracy

    Inputs:
    - y: list, ground truth
    - scores: array or list, of indices
    - map_dict: dictionary, converts indices to celltypes
    '''
    yhat = [map_dict[i] for i in scores]
    return accuracy_score(y, yhat)
