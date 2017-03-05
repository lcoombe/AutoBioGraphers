#PathMatch code
import networkx as nx

def read_path(path_filename, path):
    with open(path_filename, 'r') as f:
        for line in f:
            path.append(line.rstrip())

def read_network(network_filename, G):
    with open(network_filename, 'r') as f:
        for line in f:
            splits = line.rstrip().split()
            if len(splits) > 1:
                x = Node(xxx)
                G.add_node(x)
                #The line v1 v2 v3 v4 defines edges (v1,v2), (v1,v3) and (v1,v4).
                for n in splits[1:]
                    y = Node(xxx)
                    G.add_node(y)
                    G.add_edge(x, y)

    ## TODO CHECK IF GRAPH DIRECTED BY SEEING IF edge(x,y) and edge(y, x) exist


# Read in the correspondence information 
def read_corr(corrfilename, Gp):
    with open(corr_filename, 'r') as f:
        for line in corr_filename:
            splits = line.rstrip().splits()
            assert(len(splits)==3)
            # Pair = v1, v2, id1, id2, sim score
            p = Pair(splits[0], splits[1], splits[2])
            Gp.add_node(p)

def calculate_wieghts(Gp, G, path):
    #For each pair of vertices vij and vi+d,l satisfying 0 < d ≤ m + 1, compute the length of the shortest path d′ from vij to vi+d,l in G, and construct a directed edge from vij to vi+d,l if d′ ≤ m + 1. Add a directed edge from s to each vij and from each vij to t.

    #Each edge (vij , vi+d,l ) represents that it is always possible to find a path from vij to vi+d,l in G so that there are a total of at most m mismatches or indels between the matches (pi,vij) and (pi+d,vi+d,l). To impose mismatch and indel penalties, if ∆ ≤ 0 represents both the mismatch and indel penalty, set the weight of each edge (vij , vi+d,l ) to (max(d, d′ ) − 1)∆, the weight of each edge (s, vij ) to (i − 1)∆, and the weight of each edge (vij , t) to (n − i)∆. The above construction reduces the path matching problem to finding a path p′ from s to t in G′ with the maximum sum of vertex and edge weights.

    for x in Gp.nodes():
        # 0 based indexing p = p0, p1, p2
        i = path.index(x.v1)
        for y in Gp.nodes():
            i_d = path.index(y.v1)
            d = i_d - i
            if 0 < d <= m+1:
                dp = nx.shortest_path(G, x.v2, y.v2)
                if dp <= m+1:
                    w = (max(dp, d)-1)*delta
                    Gp.add_edge(x, y, weight=w)

        # Start
        start_weight = (i-1)*delta
        Gp.add_edge("start", x, weight=start_weight)

        # End
        end_weight = (Gp.number_of_nodes()-i)*delta
        Gp.add_edge("end", x, weight=end_weight)

def find_k_highest_scoring_paths(Gp):
    #When we assume that each edge in G′ represents only mismatches and indels and ignore different variations of mismatches or indels that can appear in a path, we can find a set of k highest scoring paths in G′ by reducing the problem to finding k shortest paths from s to t in the following modified graph with edge weights only: first negate all the vertex and edge weights in G′, then move each vertex weight into all its outgoing edges by changing each w(u, v) to w(u) + w(u, v) and setting w(v) = 0 for all vertices v. Note that the modified graph may have negative edge weights.

    # Construct Gpp 
    # find k shortest paths
    


def get_args():
    parser = argparse.ArgumentParser(description='Pathmatch')
    parser.add_argument('q', type=str, help='file name containing query path')
    parser.add_argument('i', type=str, help='file name containing input graph')
    parser.add_argument('c', type=str, help='file name containing correspondences between vertices')
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    
    # Read query path
    path = []
    read_path(args.q, path)

    # Read and create graph G from network file
    G = nx.Graph()
    read_network(network_filename, G)

    # Read and create graph G' from correspondence file
    Gp = nx.Graph()
    read_corr(corr_filename, Gp)

    # Calculate the edge weights in G'
    calculate_weights(Gp, G, p)

    # Find k highest scoring paths in G
    find_k_highest_scoring_paths(Gp)


