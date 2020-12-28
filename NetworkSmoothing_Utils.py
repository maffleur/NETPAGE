#!/usr/bin/python
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import random
import matplotlib.pyplot as plt


def MapToHGNC(net, scope):
    """
    Given a graph pandas DataFrame in the form of the networks available at HumanBase
    (http://hb.flatironinstitute.org/download) or STRING (https://string-db.org/cgi/download.pl?sessionId=Qg0Q4Hoo49KH&species_text=Homo+sapiens),
    convert the gene IDs in the first two columns to HGNC symbols
    and return a network look-up table.

    :param net: pandas graph DataFrame [gene1, gene2, Edge]
    :param scope: query scope, i.e., naming convention used for genes in the network file - to be mapped to HGNC;
                  accepted values are 'entrezgene', 'ensembl.gene', or 'ensembl.protein'
    :return: net_LUT [gene1, gene2, Edge, HGNC gene1, HGNC gene2]
    """
    # These conversions are necessary to ensure that the merging below works
    net[0] = net[0].apply(str)
    net[1] = net[1].apply(str)

    # Check the input DataFrame structure
    if net.shape[1] != 3:
        sys.exit('Error: input graph DataFrame has the wrong number of columns')

    # If the test is passed:
    else:
        import mygene
        mg = mygene.MyGeneInfo()

        NetworkGeneSymbols1 = mg.querymany(np.unique(net[:][0]), scopes=scope, fields='symbol', species='human',
                                       as_dataframe=True)
        NetworkGeneSymbols1['query'] = NetworkGeneSymbols1.index

        NetworkGeneSymbols2 = mg.querymany(np.unique(net[:][1]), scopes=scope, fields='symbol', species='human',
                                       as_dataframe=True)
        NetworkGeneSymbols2['query'] = NetworkGeneSymbols2.index

        NetworkGeneSymbols1 = NetworkGeneSymbols1.fillna(value=-9)
        NetworkGeneSymbols2 = NetworkGeneSymbols2.fillna(value=-9)

        net_LUT = pd.merge(left=net,
                           right=NetworkGeneSymbols1[['query', 'symbol']],
                           how='left',
                           left_on=net[0],
                           right_on=NetworkGeneSymbols1['query'])

        net_LUT = pd.merge(left=net_LUT,
                           right=NetworkGeneSymbols2[['query', 'symbol']],
                           how='left',
                           left_on=net_LUT[1],
                           right_on=NetworkGeneSymbols2['query'])

        net_LUT = net_LUT.iloc[:, [0, 1, 2, 4, 6]]
        net_LUT = net_LUT.rename(columns={0: "Gene1",
                                              1: "Gene2",
                                              2: "Edge",
                                              "symbol_x": "GeneSymbol1",
                                              "symbol_y": "GeneSymbol2"})
        return net_LUT


def SelectCommonGenes(network_LUT, geneburden):
    """
    Subsets the network look-up table and the gene burden table to the same set of genes in common

    :param network_LUT: pandas DataFrame [Entrez gene1, Entrez gene2, Edge, HGNC gene1, HGNC gene2]
    :param geneburden: pandas DataFrame [patients x genes]
    :return: the same DataFrames in input, but subsetted to the intersection of genes
    """
    NGSgenes = geneburden.axes[1]
    NetworkGenes = np.union1d(pd.unique(network_LUT['GeneSymbol1']), pd.unique(network_LUT['GeneSymbol2']))
    nodes = np.intersect1d(NGSgenes, NetworkGenes)

    geneburden = geneburden.loc[:, nodes]
    network_LUT = network_LUT[network_LUT['GeneSymbol1'].isin(nodes) & network_LUT['GeneSymbol2'].isin(nodes)]

    return network_LUT, geneburden, nodes


def BinariseNetwork(lut, N):
    """
    Sorts the edge values in the network LUT in descending order;
    Then assigns a 1 to the top N% edges, and 0 to the rest

    :param lut: network LUT
    :param N: percentage of top edges desired
    :return: lut with an added indicator variable [1 for top edges, 0 for remaining]
    """
    lut_sorted = lut.sort_values(by='Edge', ascending=False)

    topN = np.int(np.ceil(lut_sorted.shape[0]*N/100))
    lastN = lut_sorted.shape[0] - topN

    lut_sorted['indicator'] = np.concatenate((np.ones(topN), np.zeros(lastN))).astype(int)
    # lut_sorted.assign(indicator=np.concatenate((np.ones(topN), np.zeros(lastN))).astype(int))
    lut = lut_sorted.sort_index()

    return lut


def RestrictToPathwayGenes(net, nodes, genelist):
    """
    
    :param net: numpy matrix (adjacency!) of graph
    :param nodes: node HGNC names (same order as net rows)
    :param genelist: list of genes in the pathway (HGNC)
    :return: 
    """
    import networkx as nx
    G = nx.from_numpy_matrix(net)  # This works only with an adjacency matrix

    mapping = dict(zip(np.arange(len(nodes)), nodes))
    G = nx.relabel_nodes(G=G, mapping=mapping)

    g = nx.empty_graph()
    for gene in genelist:
        if gene in G.nodes():
            subg = nx.ego_graph(G=G, n=gene, radius=1)
            g = nx.compose(g, subg)

    S = nx.to_numpy_matrix(G=g)
    return np.asarray(S), list(g.nodes)



def DegreeNormaliseNetwork(net):
    """
    Normalises an input graph adjacency or weights matrix by the node degrees or node strengths, respectively

    :param net: graph adjacency or weights matrix (square ndarray)
    :return: A: input ndarray normalised by node degree (if adjacency) or by node strength (if weights)
    """
    net = np.asarray(net) # Just for extra safety
    net.astype(int)
    if net.shape[0] != net.shape[1]:
        sys.exit('Error: input should be a square adjacency or weights matrix')
    elif np.count_nonzero(net) == 0:
        sys.exit('Error: input network does not contain any valid edges')
    else:
        v = (np.sum(net, axis=0) ** -1)
        v[np.isinf(v)] = 0
        D = np.diag(v, k=0)
        # After Vanunu et al, PLoS Comp Bio 2010
        # A = np.dot(np.sqrt(D), np.dot(net, np.sqrt(D)))
        # After Hofree et al., Nat Meth 2013; faster than line above (only one dot prod), but same eigenvalues
        A = np.dot(D, net)

    return A


def quantile_norm(X):
    """
    From https://github.com/elegant-scipy/elegant-scipy/blob/master/markdown/ch2.markdown
    Normalize the columns of X to each have the same distribution.

    Given an expression matrix (microarray data, read counts, etc) of M genes
    by N samples, quantile normalization ensures all samples have the same
    spread of data (by construction).

    The data across each row are averaged to obtain an average column. Each
    column quantile is replaced with the corresponding quantile of the average
    column.

    Parameters
    ----------
    X : 2D array of float, shape (M, N)
        The input data, with M rows (genes/features) and N columns (samples).

    Returns
    -------
    Xn : 2D array of float, shape (M, N)
        The normalized data.
    """
    # compute the quantiles
    quantiles = np.mean(np.sort(X, axis=0), axis=1)

    # compute the column-wise ranks. Each observation is replaced with its
    # rank in that column: the smallest observation is replaced by 1, the
    # second-smallest by 2, ..., and the largest by M, the number of rows.
    ranks = np.apply_along_axis(stats.rankdata, 0, X)

    # convert ranks to integer indices from 0 to M-1
    rank_indices = ranks.astype(int) - 1

    # index the quantiles for each rank with the ranks matrix
    Xn = quantiles[rank_indices]

    return Xn


def TableAndTest(mutationstatus, casecontrol):
    contingency = pd.crosstab(mutationstatus, casecontrol)
    if contingency.shape != (2, 2):
        return np.nan, np.nan
    else:
        oddsratio, pvalue = stats.fisher_exact(contingency)
        return oddsratio, pvalue


def void_protectives(geneburden, outcome):

    testres = geneburden.apply(func=TableAndTest, args=[outcome], axis=0, raw=False)
    testres = testres.apply(pd.Series)
    testres.columns = ['Effect', 'Pvalue']
    #protectives = testres[(testres['Effect'] < 1) & (testres['Pvalue'] < 0.01)].index
    protectives = testres[testres['Effect'] < 1].index
    geneburden[protectives] = 0

    return geneburden


def void_risks(geneburden, outcome):

    testres = geneburden.apply(func=TableAndTest, args=[outcome], axis=0, raw=False)
    testres = testres.apply(pd.Series)
    testres.columns = ['Effect', 'Pvalue']
    #risks = testres[(testres['Effect'] > 1) & (testres['Pvalue'] < 0.01)].index
    risks = testres[testres['Effect'] > 1].index
    geneburden[risks] = 0

    return geneburden


def invert_protectives(geneburden, outcome):

    testres = geneburden.apply(func=TableAndTest, args=[outcome], axis=0, raw=False)
    testres = testres.apply(pd.Series)
    testres.columns = ['Effect', 'Pvalue']
    # protectives = testres[(testres['Effect'] < 1) & (testres['Pvalue'] < 0.01)].index
    protectives = testres[testres['Effect'] < 1].index
    geneburden[protectives] = -geneburden[protectives]

    return geneburden


def logOR(geneburden, outcome):

    testres = geneburden.apply(func=TableAndTest, args=[outcome], axis=0, raw=False)
    testres = testres.apply(pd.Series)
    testres.columns = ['Effect', 'Pvalue']
    testres['Beta'] = np.log(testres['Effect'])
    # replace with log(Effect)
    for gene in testres.index:
        geneburden[gene] = testres.loc[gene, 'Beta'] * geneburden[gene]

    return geneburden


def none(geneburden, outcome):
    return geneburden


def AnalyseDirectionOfEffects(geneburden, outcome, effect):
    """

    :param geneburden: gene columns in the merged dataframe
    :param outcome: outcome column in the merged dataframe
    :param effect: choices=['void_protectives', 'void_risks', 'invert_protectives', 'logOR', 'none']
    :return:
    """
    # Outer switch on the values of argument 'effect'
    # In each case in the switch statement, nest an if clause on the value of typ
    # In both branches of the if clause, tests are to be run column-wise in a vectorised way through pandas.DataFrame.apply

    switcher = {
        'void_protectives': void_protectives,
        'void_risks': void_risks,
        'invert_protectives': invert_protectives,
        'logOR': logOR,
        'none': none
    }
    geneburden_modified = switcher[effect](geneburden, outcome)
    return geneburden_modified


def randomize_by_edge_swaps(G, num_iterations):

    newgraph = G.copy()
    edge_list = list(newgraph.edges())
    num_edges = len(edge_list)
    total_iterations = num_edges * num_iterations

    for i in xrange(total_iterations):
        rand_index1 = int(round(random.random() * (num_edges - 1)))
        rand_index2 = int(round(random.random() * (num_edges - 1)))
        original_edge1 = edge_list[rand_index1]
        original_edge2 = edge_list[rand_index2]
        head1, tail1 = original_edge1
        head2, tail2 = original_edge2

        # Flip a coin to see if we should swap head1 and tail1 for
        # the connections
        if random.random() >= 0.5:
            head1, tail1 = tail1, head1

        # The plan now is to pair head1 with tail2, and head2 with
        # tail1
        #
        # To avoid self-loops in the graph, we have to check that,
        # by pairing head1 with tail2 (respectively, head2 with
        # tail1) that head1 and tail2 are not actually the same
        # node. For example, suppose we have the edges (a, b) and
        # (b, c) to swap.
        #
        #   b
        #  / \
        # a   c
        #
        # We would have new edges (a, c) and (b, b) if we didn't do
        # this check.

        if head1 == tail2 or head2 == tail1:
            continue

        # Trying to avoid multiple edges between same pair of nodes;
        # for example, suppose we had the following
        #
        # a   c
        # |*  |           | original edge pair being looked at
        # | * |
        # |  *|           * existing edge, not in the pair
        # b   d
        #
        # Then we might accidentally create yet another (a, d) edge.
        # Note that this also solves the case of the following,
        # missed by the first check, for edges (a, b) and (a, c)
        #
        #   a
        #  / \
        # b   c
        #
        # These edges already exist.

        if newgraph.has_edge(head1, tail2) or newgraph.has_edge(
                head2, tail1):
            continue

        # Suceeded checks, perform the swap
        original_edge1_data = newgraph[head1][tail1]
        original_edge2_data = newgraph[head2][tail2]

        newgraph.remove_edges_from((original_edge1, original_edge2))

        new_edge1 = (head1, tail2, original_edge1_data)
        new_edge2 = (head2, tail1, original_edge2_data)
        newgraph.add_edges_from((new_edge1, new_edge2))

        # Now update the entries at the indices randomly selected
        edge_list[rand_index1] = (head1, tail2)
        edge_list[rand_index2] = (head2, tail1)

    assert len(newgraph.edges()) == num_edges
    return newgraph


def network_propagate(GeneBurden, network, alpha, tol, quantile):
    """

    :param GeneBurden: self-explanatory
    :param network: self-explanatory
    :param alpha: diffusion length (only one value allowed)
    :param tol: convergence tolerance
    :param quantile: True/False
    :return:
    """
    # print 'Propagating the gene burden through the network via network smoothing...'

    F_0 = GeneBurden.values
    A = DegreeNormaliseNetwork(network)


    print('Network propagation with diffusion length ' + str(alpha))
    F_t = F_0
    diff = 1
    maxdiff = 1e5
    it = 1
    while diff > tol:
        if diff >= maxdiff:
            sys.exit("Error: network propagation did not converge")
        else:
            print 'Iteration no. ' + str(it)
            F_T = alpha * np.dot(F_t, A) + (1 - alpha) * F_0
            diff = np.linalg.norm(F_T - F_t, ord=2)
            print 'Convergence difference: ' + str(diff)
            F_t = F_T
            it += 1
    if quantile:
        q = "quantiles"
        print 'Quantile-normalising the original and smoothed profiles...'
        SmoothProfile_quantiles = quantile_norm(F_T.transpose())
        SmoothProfile_quantiles = SmoothProfile_quantiles.transpose()
    else:
        q = "NOquantiles"
        SmoothProfile_quantiles = F_T

    return SmoothProfile_quantiles, q


def remove_edges(G, p_remove_connection):
    '''
    for each node,
    remove a connection, with prob p_remove_connection

    operates on G in-place

    Adapted from https://stackoverflow.com/questions/42591549/add-and-delete-a-random-edge-in-networkx
    '''

    rem_edges = []

    for node in G.nodes():
        # find the other nodes this one is connected to
        connected = [to for (fr, to) in G.edges(node)]

        # probabilistically remove a random edge
        if len(connected): # only try if an edge exists to remove
            if random.random() < p_remove_connection:
                remove = random.choice(connected)
                G.remove_edge(node, remove)
                # print "\tedge removed:\t {} -- {}".format(node, remove)
                rem_edges.append( (node, remove) )
                # book-keeping, in case lists are important later?
                connected.remove(remove)
    return rem_edges
