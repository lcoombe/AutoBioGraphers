# coding: utf-8
#PathMatch code

import argparse
import math
import networkx as nx

from ep import find_k_shortest_paths

INDEL_PENALTY = -1
MAX_NO_OF_MISMATCHES_AND_GAPS = 1
NO_OF_TOP_PATHS_OUTPUT = 3
CORR_E_VALUE_CUTOFF = float(1e200)

def safe_log(sim):
    if sim < 0:
        print "Negative sim score"
        sys.exit(1)
    elif sim <= 1e-299:
        return math.log(1e-299)
    else:
        return math.log(sim)

# Read query path from file and load into list path
def read_path(path_filename, path):
    with open(path_filename, 'r') as f:
        for line in f:
            path.append(line.rstrip())

# Read input graph from file and construct input graph G
def read_network(network_filename, G):
    with open(network_filename, 'r') as f:
        for line in f:
            #The line v1 v2 v3 v4 defines edges (v1,v2), (v1,v3) and (v1,v4).
            splits = line.rstrip().split()
            if len(splits) > 1:
                x = splits[0]
                G.add_node(x)
                for y in splits[1:]:
                    G.add_node(y)
                    G.add_edge(x, y)


# Read correspondence information from file and construct graph G' nodes only including v1 in query path p
def read_corr(corr_filename, path, Gp):
    with open(corr_filename, 'r') as f:
        for line in f:
            # [v1 v2 corr] where v1 is a vertex in the query path, v2 is 
            # a vertex in the input graph and corr is the similarity score
            splits = line.rstrip().split()
            assert(len(splits)==3)

            if splits[0] not in path:
                continue

            # p = (pi, vij)
            p = (splits[0], splits[1])
            # -l= "type of match score"
            #       1 -- use negative logarithm of similarity score as match score
            w = -safe_log(float(splits[2]))
            if float(splits[2]) <= CORR_E_VALUE_CUTOFF:
                Gp.add_node(p, weight=w)

# Calculate the edge weights of graph G'
def calculate_weights(Gp, G, path):
    # From the paper:
    # For each pair of vertices vij and vi+d,l satisfying 0  < d <= m + 1, compute the length of the shortest path d' from vij to vi+d,l in G, and construct a directed edge from vij to vi+d,l if d' <= m + 1. Add a directed edge from s to each vij and from each vij to t.
    # Each edge (vij , vi+d,l) represents that it is always possible to find a path from vij to vi+d,l in G so that there are a total of at most m mismatches or indels between the matches (pi,vij) and (pi+d,vi+d,l). To impose mismatch and indel penalties, if DELTA <= 0 represents both the mismatch and indel penalty, set the weight of each edge (vij , vi+d,l ) to (max(d, d' ) − 1)DELTA, the weight of each edge (s, vij ) to (i − 1)DELTA, and the weight of each edge (vij , t) to (n − i)DELTA. The above construction reduces the path matching problem to finding a path p' from s to t in G' with the maximum sum of vertex and edge weights.

    # Adding source and sink nodes with weight 0
    Gp.add_node("s", weight=0)
    Gp.add_node("t", weight=0)

    # Loop through all possible pairs of verticies except for start, end nodes
    # node = (pi-vertex in query path, vij-vertex in input graph)
    for x in Gp.nodes():
        if x == "s" or x == "t":
            continue
        i = path.index(x[0]) # Find out value of i for node x = vi,j

        for y in Gp.nodes():
            if y == "s" or y == "t" or y == x:
                continue
            i_d = path.index(y[0]) # Find out value of i+d for node y = vi+d,l
            d = i_d - i
            if 0 < d <= MAX_NO_OF_MISMATCHES_AND_GAPS+1:
                if x[1] in G.nodes() and y[1] in G.nodes():
                    try:
                        # Find the shortest path d' in G from vi,j to vi+d,l
                        shortest_path = nx.shortest_path(G, x[1], y[1])
                        dp = len(shortest_path)-1 # Since source and target node are included in path
                    except Exception: # No path exists between source and target
                        dp = float("inf")

                    if 0 < dp <= MAX_NO_OF_MISMATCHES_AND_GAPS+1:
                        # Set edge weight to impose mismatch and indel penalties
                        w = (max(dp, d)-1)*INDEL_PENALTY
                        Gp.add_edge(x, y, weight=w)

        # Start
        start_weight = (i)*INDEL_PENALTY
        Gp.add_edge("s", x, weight=start_weight)

        # End
        end_weight = (len(path)-i-1)*INDEL_PENALTY
        Gp.add_edge(x, "t",  weight=end_weight)

def find_k_highest_scoring_paths(Gp):

    # From paper:
    # When we assume that each edge in G' represents only mismatches and indels and ignore different variations of mismatches or indels that can appear in a path, we can find a set of k highest scoring paths in G' by reducing the problem to finding k shortest paths from s to t in the following modified graph with edge weights only: first negate all the vertex and edge weights in G', then move each vertex weight into all its outgoing edges by changing each w(u, v) to w(u) + w(u, v) and setting w(v) = 0 for all vertices v. Note that the modified graph may have negative edge weights.

    # Construct Gpp, the modified G' graph with negated edge+node weights
    Gpp = nx.DiGraph()

    # Add all nodes from Gp to Gpp, with weight w(v) = 0
    for n in Gp.nodes(): # node = (pi-vertex in query path, vij-vertex in input graph)
        Gpp.add_node(n)

    # Adding all edges from Gp to Gpp with weight -w(u) + -w(u, v), since all node+edge weights were negated
    node_weights = nx.get_node_attributes(Gp, 'weight')
    for n in Gp.nodes():
        edges = Gp.in_edges(n, True) # Get all out edges for the node
        negated_node_weight = -node_weights[n]
        for e in edges: # e = (node, neighbor, data)
            negated_edge_weight = -e[2]['weight']
            w = negated_node_weight + negated_edge_weight 
            Gpp.add_edge(e[0], e[1], weight=w)

    return find_k_shortest_paths(Gpp, NO_OF_TOP_PATHS_OUTPUT)


def shorten_result(col1, col2, col3):

    def _helper(col1, col2, col3):
        for i in range(len(col1))[:-1]:
            a1, b1 = col1[i], col3[i]
            a2, b2 = col1[i+1], col3[i+1]

            if a1 == "-" and b1 != "-":
                if a2 != "-" and b2 == "-":
                    c1, c2 = a2, b1
                    new_col1 = col1[:i] + [c1] + col1[i+2:]
                    new_col3 = col3[:i] + [c2] + col3[i+2:]
                    new_col2 = col2[:i] + col2[i+1:]
                    return [True, new_col1, new_col2, new_col3]

            if a1 != "-" and b1 == "-":
                if a2 == "-" and b2 != "-":
                    new_col1 = col1[:i] + [c1] + col1[i+2:]
                    new_col3 = col3[:i] + [c2] + col3[i+2:]
                    new_col2 = col2[:i] + col2[i+1:]
                    return [True, new_col1, new_col2, new_col3]
        return [False, [], [], []]


    new_col1, new_col2, new_col3 = col1, col2, col3
    result = _helper(col1, col2, col3)

    while result[0]:
        flag, new_col1, new_col2, new_col3 = result
        result  = _helper(new_col1, new_col2, new_col3)
    return new_col1, new_col2, new_col3

# Format results and write to outfile given k shortest paths
# Results = [ (score, [('s', n1), (n1, n2), (n2, 't')])....]
def format_result(path, G, results, outfile):
    for index, item in enumerate(results):
        score, long_result = item
        result = [x[0] for x in long_result[1:]]
        col1, col2, col3 = [], [], [] # Each col in outfile respectively 
        result_q = [x[0] for x in result]

        for q in path:
            col1.append(q)
            if q not in result_q: # Deletion
                col2.append("")
                col3.append("-")
            else: # Hit/match
                col2.append("--")
                entry = result[result_q.index(q)]
                assert(q == entry[0])
                col3.append(entry[1])

                # Insertions
                if result_q.index(q)+1 < len(result):
                    next_entry = result[result_q.index(q)+1]
                    # Find the shortest path d' in G from current entry to next entry
                    shortest_path = nx.shortest_path(G, entry[1], next_entry[1])
                    for a in shortest_path[1:-1]:
                        col1.append("-")
                        col2.append("")
                        col3.append(a)

        assert(len(col1)+len(col2)+len(col3) == len(col1)*3)
        new_col1, new_col2, new_col3 = shorten_result(col1, col2, col3)
        assert(len(new_col1)+len(new_col2)+len(new_col3) == len(new_col1)*3)

        if index == 0:
            typ = 'w'
        else:
            typ = 'a'

        # write results to outfile
        with open(outfile, typ) as f:
            f.write("Result " + str(index+1) + ": score=" + str(score) + "\n")
            rows = [[new_col1[i], new_col2[i], new_col3[i]] for i in range(len(new_col1))]

            widths = [max(map(len, col)) for col in zip(*rows)]
            for row in rows:
                line = "    ".join((val.ljust(width) for val, width in zip(row, widths)))
                f.write(line+"\n")

            f.write("\n\n")
    
# Parse input arguments
def get_args():

    parser = argparse.ArgumentParser(description='Pathmatch')
    parser.add_argument('query', type=str, help='file name containing query path')
    parser.add_argument('input', type=str, help='file name containing input graph')
    parser.add_argument('corr', type=str, help='file name containing correspondences between vertices')
    parser.add_argument("output", type=str, help="output file name")
    parser.add_argument("-g", type=int, help="maximum number of mismatches or indels between two matches, default=1")
    #parser.add_argument("-l", type=int, choices=[0, 1], help="type of match score: 0 -- use similarity score as match score, 1 -- use negative logarithm of similarity score as match score, default=1")
    #parser.add_argument("-s", type=int, choices=[0, 1], help="whether gaps are allowed at the start and the end of path alignment: 0 -- not allowed, i.e. must be matches at the start and end, 1 -- allowed, default=1")
    parser.add_argument("-p", type=float, help="mismatch and indel penalty (>0), default=1")
    parser.add_argument("-n", type=int, help="number of output paths, default=3")

    args = parser.parse_args()
    return args

def main():
    args = get_args()

    if args.g is not None:
        global MAX_NO_OF_MISMATCHES_AND_GAPS 
        MAX_NO_OF_MISMATCHES_AND_GAPS = args.g
    if args.p is not None:
        global INDEL_PENALTY
        INDEL_PENALTY = -args.p
    if args.n is not None:
        global NO_OF_TOP_PATHS_OUTPUT 
        NO_OF_TOP_PATHS_OUTPUT = args.n
    
    # Read query path
    path = []
    read_path(args.query, path)

    # Read and create graph G from input network file
    G = nx.Graph()
    read_network(args.input, G)

    # Read and create graph G' from correspondence file
    Gp = nx.DiGraph() # directed graph
    read_corr(args.corr, path, Gp)

    # Calculate the edge weights in G'
    calculate_weights(Gp, G, path)

    # Find k highest scoring paths in G
    results = find_k_highest_scoring_paths(Gp)

    # Format and write reslts to outfile
    format_result(path, G, results, args.output)

if __name__ == "__main__":
    main()
