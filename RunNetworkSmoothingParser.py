from RunNetworkSmoothingMain import main
import argparse
import sys

parser = argparse.ArgumentParser(add_help=True, version='1.0', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i', '--inputnet',
                    action='store',
                    help='Path and name of the input network file (required)\n\n',
                    required=True)
parser.add_argument('-g', '--geneburden',
                    action='store',
                    help='Path and name of the gene burden table file (required)\n\n',
                    required=True)
parser.add_argument('--typ',
                    help='Type of rare variant encoding (options allowed: binary, burden; required)\n\n',
                    choices=['binary', 'burden'],
                    action='store',
                    required=True)
parser.add_argument('--scope',
                    help='Gene naming convention used in the network file.\n'
                         'Needed to map gene IDs to HGNC symbols.\n\n',
                    choices=['entrezgene', 'ensembl.gene', 'ensembl.protein'],
                    action='store',
                    required=True)
parser.add_argument('-o', '--outcome',
                    action='store',
                    help='File storing the outcome wrt which direction of effect must be considered (required)\n\n')
parser.add_argument('-e', '--effect',
                    help='How do you want to deal with direction of effects? (case-sensitive, default "none")\n'
                         '* "void_protectives": set to 0 the mutation status/burden of protective genes\n'
                         '* "void_risks": set to 0 the mutation status/burden of risk genes\n'
                         '* "inversion": set mutation status to -1 for protective, +1 for risk genes\n'
                         '* "logOR": replace mutation status with logOR from Fisher test\n'
                         '* "none": do nothing\n\n',
                    default='none',
                    action='store',
                    choices=['void_protectives', 'void_risks', 'invert_protectives', 'logOR', 'none'])
parser.add_argument('-t', '--tol',
                    action='store',
                    type=float,
                    help='Convergence threshold for network smoothing (float, default 1e-6)\n\n',
                    default=1e-6)
parser.add_argument('-a', '--alpha',
                    # nargs='*',
                    help='Diffusion distance allowed for a mutation signal through the network (float, default 0.5).\n'
                         'Multiple values not allowed\n\n',
                    default=0.5,
                    type=float)
parser.add_argument('-s', '--selfloops',
                    action='store_true',
                    help="Should the network have self loops? (default False)\n\n")
parser.add_argument('--binarynet',
                    action='store_true',
                    help='Binarize the weighted network to an adjacency matrix;\n'
                         'pick only the top p%% edges, as supplied with -p (default False)\n\n')
parser.add_argument('-p', '--percentage',
                    action='store',
                    help='Percentage of top edges to retain when binarising the network (float, default 1)\n\n',
                    type=float,
                    default=1,
                    required='--binarynet' in sys.argv)
parser.add_argument('--edgeremoval',
                    action='store_true',
                    help='For each gene, remove a random edge with probability p-rem,'
                         'as supplied with --prob_removal (default False)\n\n')
parser.add_argument('--probremoval',
                    action='store',
                    help='Probability of randomly removing an edge incident on a gene (default 0.1)\n\n',
                    type=float,
                    default=0.1,
                    required='--edgeremoval' in sys.argv)
parser.add_argument('-r', '--randomnet',
                    action='store_true',
                    help="If present, perform degree-preserving randomisation of the network after binarisation\n\n")
parser.add_argument('--nrand', '-n',
                    action='store',
                    type=int,
                    default=30,
                    help='How many replicates of the randomised network? (default 30)\n\n',
                    required='--randomnet' in sys.argv)
parser.add_argument('--quantile', '-q',
                    action='store_true',
                    help="Quantile-normalise the smooth profile's rows (default False)\n\n")
parser.add_argument('--pathway', '-pw',
                    action='store',
                    help='Restrict network and mutation burden to genes in a pathway.\n'
                         'Please provide name of the file containing the list of genes in the desired pathway\n'
                         '(one gene per line)\n\n')
parser.add_argument('--plots', '-pl',
                    action='store_true',
                    help='Plot a heatmap of the network? (default False)\n\n')

args = parser.parse_args()


main(args.inputnet, args.geneburden, args.typ, args.scope, args.outcome, args.effect, args.tol, args.alpha,
     args.selfloops, args.binarynet, args.percentage, args.edgeremoval, args.probremoval,
     args.randomnet, args.nrand,
     args.quantile, args.pathway, args.plots)
