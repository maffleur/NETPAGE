from NetworkSmoothing_Utils import *
import os
import networkx as nx
import matplotlib.pyplot as plt
import sys


def main(f1, f2, typ, scope, f3, effect='invert_protectives', tol=1e-6, alpha=0.5, selfloops=False, binarynet=True, percentage=1, edgeremoval=False, probremoval=0.1, randomnet=True, nrand=30, quantile=False, pathway="", plots=False):

    dir_path = os.path.dirname(os.path.realpath(f2))
    gene_base = os.path.basename(os.path.realpath(f2))
    net_base = os.path.basename(os.path.realpath(f1))

    print 'Reading network file...'
    network = pd.read_table(f1, header=None)
    if network.shape[1] == 2 & binarynet:
        sys.exit('Error: You are trying to binarise an unweighted network. This is not allowed. '
                 'Please remove the --binarynet flag and try again.')
    if network.shape[1] == 2:
        network.insert(2, 'Edge', 1)

    print 'Reading gene burden table file...'
    GeneBurden = pd.read_table(f2)
    GeneBurden['PTID'] = GeneBurden['PTID'].str.upper()
    subjects = GeneBurden['PTID']

    if typ == 'binary':
        print 'Reading outcome (phenotype) file...'
        outcome = pd.read_table(f3)
        print 'Merging gene burden and outcome tables...'
        GeneBurden = pd.merge(left=outcome,
                              right=GeneBurden,
                              how='right',
                              on='PTID')
        GeneBurden = GeneBurden.drop_duplicates()
        print 'Dealing with direction of mutation effects...'
        GeneBurden = AnalyseDirectionOfEffects(GeneBurden.iloc[:, 2:], GeneBurden.iloc[:, 1], effect)

    print 'Converting network gene IDs to HGNC symbols...'
    # Use scope='entrezgene' with GIANT networks,
    # scope='ensembl.protein' with STRING network,
    # scope='ensembl.gene' with HuRI network
    network_LUT = MapToHGNC(network, scope=scope)

    print 'Selecting common genes...'
    network_LUT, GeneBurden, nodes = SelectCommonGenes(network_LUT, GeneBurden)

    # If required, binarize the network by choosing the top p% edges
    if binarynet:
        b = 'binary'
        print 'Picking the top ' + str(percentage) + '% edges for the adjacency matrix...'
        network_LUT_bin = BinariseNetwork(network_LUT, percentage)
        print(network_LUT_bin.head(n=20))

        print 'Reshaping matrix...'
        squarenet = network_LUT_bin.pivot_table(index='GeneSymbol1', columns='GeneSymbol2', values='indicator')
        squarenet = squarenet.reindex(index=nodes, columns=nodes, fill_value=0)

    else:
        if max(network_LUT['Edge']) == 1:
            b = 'binary'  # make it work for HuRI, which is already binary so doesn't go through the if branch above
        else:
            b = 'weighted'
        # r = "NotRandomised"  # This assumes that randomisation is not allowed on a weighted graph
        print 'Reshaping matrix...'
        squarenet = network_LUT.pivot_table(index='GeneSymbol1', columns='GeneSymbol2', values='Edge')
        squarenet = squarenet.reindex(index=nodes, columns=nodes, fill_value=0)

    # Use an ndarray to make the network symmetric and to fill the network diagonal with 0s or 1s
    N = squarenet.values
    N[np.isnan(N)] = 0  # Necessary to get rid of residual NaNs
    N = np.triu(N) + np.tril(N).T + np.tril(N) + np.triu(N).T

    if selfloops:  # Should diagonal be filled with zeros or ones?
        np.fill_diagonal(N, 1)
    else:
        np.fill_diagonal(N, 0)

    if randomnet:
        r = "randomised"
        print 'Generating {0} degree-preserving randomised versions of the network...'.format(nrand)
        # convert numpy array to networkx graph object
        tmp_graph = nx.from_numpy_matrix(N)
        N = []
        for n in np.arange(nrand):
            # apply randomize_by_edge_swaps() with 10 iterations
            tmp_graph_random = randomize_by_edge_swaps(tmp_graph, num_iterations=10)
            # reconvert graph object to numpy matrix
            N.append(nx.to_numpy_matrix(tmp_graph_random))
    else:
        r = "NotRandomised"

    # print 'Writing network matrix to csv...'
    # N2 = pd.DataFrame(N, index=nodes, columns=nodes)
    # f_name = os.path.join(dir_path, '{0}_{1}_{2}_NetworkMatrix_{3}.csv'.format(net_base, b, percentage, r))
    # N2.to_csv(f_name, na_rep="NA")

    if pathway:
        if randomnet:
            sys.exit('Error: does not make sense to restrict a randomised network to the pathway.'
                     'Please remove either the --randomnet or the --pathway flag and try again.\n\n')
        else:
            print 'Computing the subgraph induced by genes in the supplied pathway and their first neighbours...'
            genelist = open(pathway, 'r').read().splitlines()
            N, nodes = RestrictToPathwayGenes(net=N, nodes=nodes, genelist=genelist)
            GeneBurden = GeneBurden.loc[:, nodes]
            p = 'pathway'
    else:
        p = 'NOpathway'

    while True:
        if plots and randomnet:
            print('Will not plot all the {0} randomised versions of the network. '
                  'Sorry, but it\'s for the best.'.format(nrand))
            break
        elif plots:
            # Plot the network matrix as heatmap
            print 'Plotting the network for inspection...'
            fig, ax = plt.subplots()
            cax = ax.imshow(N, cmap='Reds', interpolation='None')
            fig.colorbar(cax, ticks=[0, 0.5, 1])

            f_name = os.path.join(dir_path, '{0}_{1}_{2}_NetworkMatrix_{3}_{4}.png'.format(net_base, b, percentage, p, r))
            plt.savefig(f_name, dpi=300)
            break

        break

    if edgeremoval and randomnet:
        sys.exit('Error: --edgeremoval is not meant to be used with --randomnet.'
                 'Please remove one of the two flags and relaunch.\n\n')
    elif edgeremoval:
        er = 'EdgesRemoved_prob{}'.format(probremoval)
        print 'Removing an edge for each node with probability {}...'.format(probremoval)
        tmp_graph = nx.from_numpy_matrix(N)
        rem_edges = remove_edges(tmp_graph, probremoval)
        N = nx.to_numpy_matrix(tmp_graph)
    else:
        er = 'EdgesNotRemoved'



    # N is the numpy array representation of the squarenet pandas df,
    # in a square SYMMETRIC format [e.g., nodes x nodes];
    # edges are the posterior probabilities for an interaction;
    # unknown edges are set to 0 and diagonal elements are set to 0 or 1 depending on user input;
    # GeneBurden is the gene burden in a rectangular DataFrame format [e.g., patients x nodes]
    # And the nodes are the same in both DataFrames

    if randomnet:
        for n in np.arange(nrand):
            print 'Randomised network no. ' + str(n+1)
            SmoothProfile_quantiles, q = network_propagate(GeneBurden=GeneBurden, network=N[n], alpha=alpha,
                                                           tol=tol, quantile=quantile)
            SmoothedGeneBurden = pd.DataFrame(SmoothProfile_quantiles, index=subjects, columns=GeneBurden.axes[1])
            f_name = os.path.join(dir_path,
                                  '{0}_{1}_{2}_{3}_{4}_{5}Repl{6}_SmoothedProfiles_{7}_alpha{8}_{9}.csv'.format(gene_base,
                                                                                                         net_base, b,
                                                                                                         percentage, p,
                                                                                                         r, n+1, q,
                                                                                                         alpha, effect))
            SmoothedGeneBurden.to_csv(f_name, na_rep="NA")
    else:
        SmoothProfile_quantiles, q = network_propagate(GeneBurden=GeneBurden, network=N, alpha=alpha,
                                                       tol=tol, quantile=quantile)
        SmoothedGeneBurden = pd.DataFrame(SmoothProfile_quantiles, index=subjects, columns=GeneBurden.axes[1])
        f_name = os.path.join(dir_path,
                              '{0}_{1}_{2}_{3}_{4}_{5}_{6}_SmoothedProfiles_{7}_alpha{8}_{9}.csv'.format(gene_base,
                                                                                                     net_base, b,
                                                                                                     percentage, p, r,
                                                                                                     er, q, alpha,
                                                                                                     effect))
        SmoothedGeneBurden.to_csv(f_name, na_rep="NA")


if __name__ == '__main__':
    f1 = 'GIANTnetworks/hippocampus_top'
    f2 = 'ADSP_WES/ADSP_GeneBurden_CADD20_binary_noADNI.txt'
    f3 = 'ADSP_WES/ADSP_noADNI_CaseControlStatus.txt'

    main(f1=f1, f2=f2, typ='binary', scope='entrezgene', f3=f3, effect='invert_protectives', alpha=0.5, binarynet=True,
         percentage=1, randomnet=True, nrand=1, edgeremoval=False, plots=False)
