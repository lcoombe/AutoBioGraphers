#GraphMatch Code

import networkx as nx
import argparse
import re
from Corr_Vertex import Corr_Vertex
from Tree import Tree
from networkx.drawing.nx_pydot import write_dot
import copy
from math import log
from Alignment import Alignment

#EXAMPLE USAGE: python GraphMatch.py test-data/query test-data/input test-data/corr_short

#GLOBAL VARIABLES
non_assoc_vertex_penalty = 1.0
indel_penalty = 0.1
m = 1
T = Tree()
results = []
scores = []
score_method = 1
max_score = 0
num_query_nodes = 0
top_alignments = []

#HELPER FUNCTIONS

#Reading in input and query graphs
def readGraph(filename):

    G_Graph = nx.Graph()
    refFile = open(filename, 'r')
    for line in refFile:
        line = line.strip()
        split_line = re.split(r'\s+', line)
        for node in split_line:
            if not G_Graph.has_node(node):
                G_Graph.add_node(node)
        for i in range(1, len(split_line)):
            if not G_Graph.has_edge(split_line[0], split_line[i]):
                G_Graph.add_edge(split_line[0], split_line[i])

    refFile.close()
    return G_Graph

# Load in correspondances between vertices in G and G0
#Format of input file: v1 v2 corr: where v1 is a vertex in the query graph, v2 is a vertex in the input graph
def readCorrespondances(corrFileName):
    corr = {} #Dictionary: key = query node name. Value = list of corresponding vertices
    corrFile = open(corrFileName, 'r')
    for line in corrFile:
        line = line.strip()
        split_line = re.split(r'\s+', line)
        query_vertex = split_line[0]
        ref_vertex = split_line[1]
        score = split_line[2]
        corresp_vertex = Corr_Vertex(ref_vertex, score, query_vertex)
        if query_vertex not in corr:
            corr[query_vertex] = []
        corr[query_vertex].append(corresp_vertex)
    corrFile.close()
    return corr

# Making G' --> will be used when enumerating solution graphs G0'
def makeGraph_prime(corr, G0, G):
    V_prime = []
    G_prime = nx.Graph()
    #Taking V' --> union of all corresponding vertices
    for queryNode in corr:
        V_prime.extend(corr[queryNode])
    #Building up the G' graph
    for i in range(0, len(V_prime)):
        vij = V_prime[i]
        for j in range(i+1, len(V_prime)):
            vkl = V_prime[j]
            if vij.queryVertex != vkl.queryVertex: #Only would consider adding edge to G' if they have different corresponding vertices
                has_g0_edge = G0.has_edge(vij.queryVertex, vkl.queryVertex) #Is there a (vi, vk) edge in G0?
                path_length_G = float("inf")
                if nx.has_path(G, vij.name, vkl.name):
                    short_path_G = nx.shortest_path(G, vij.name, vkl.name) #Path length from vij to vkl in G
                    path_length_G = len(short_path_G)-1
                if has_g0_edge and path_length_G <= m+1:
                    vij_string = vij.stringifyVertex()
                    vkl_string = vkl.stringifyVertex()
                    G_prime.add_edge(vij_string, vkl_string)
    return G_prime

#Given a vertex vi, and a list of corresponding vertex, return the corresponding vertex of vi
def getCorrespondingVertex(v, correspList):
    for c in correspList:
        if c.queryVertex == v:
            return c
    return None #No corresponding vertex

#Scoring the alignment between two graphs G0 and G0'
def ScoreAlignment(W_prime, G, G0):
    score = 0
    assoc_vertices = []
    indel_vertices = []
    query_edges = []
    ref_edges = []
    if score_method: #Use -log(corresp score) as match value. want to maximize this score.
        for node in G0.node:
            corresp = getCorrespondingVertex(node, W_prime)
            if corresp is None:
                score = score - non_assoc_vertex_penalty
            else:
                assoc_vertices.append(node)
                if corresp.score <= 1e-200:
                    score += -1.0*log(1e-200)
                else:
                    score += -1.0*log(corresp.score)

    else: #Using raw correspondance score for matching -- want to maximize this score.
        for node in G0.node:
            corresp = getCorrespondingVertex(node, W_prime)
            if corresp is None:
                score = score - non_assoc_vertex_penalty
            else:
                assoc_vertices.append(node)
                score -= corresp.score

    #Indel penalties
    for i in range(0, len(assoc_vertices)):
        for j in range(i+1, len(assoc_vertices)):
            vi = assoc_vertices[i]
            vk = assoc_vertices[j]
            if G0.has_edge(vi, vk):
                query_edges.append((vi, vk))
                vij = getCorrespondingVertex(vi, W_prime)
                vkl = getCorrespondingVertex(vk, W_prime)
                path_G = nx.shortest_path(G, vij.name, vkl.name)
                pathLen = len(path_G) - 2
                for h in range(0, pathLen):
                    indel_vertices.append(path_G[h+1])
                score = score - pathLen*indel_penalty
                ref_edges.append(path_G)

    return Alignment(score, copy.copy(W_prime), indel_vertices, query_edges, ref_edges)


#Determining if the vertex combos are a valid solution
def isValidSolution(V0_plus, W_prime, G0, G_prime):
    for vi in V0_plus:
        for vk in V0_plus:
            if vi != vk and G0.has_edge(vi, vk):
                vij = getCorrespondingVertex(vi, W_prime)
                vkl = getCorrespondingVertex(vk, W_prime)
                if not G_prime.has_edge(vij.stringifyVertex(), vkl.stringifyVertex()):
                    return False
    return True

#Returns True if based on the branch and bound approach, the recurrence should continue. Returns False if can't get to better alignment than top k
# alignments down this branch, so the search is pruned.
def continueDFS(W, score, lower_bound):
    num = len(W)
    if score_method:
        l = -1*log(1e-200)
    else:
        return True
    upper_bound = score + l*(num_query_nodes - num)
    if upper_bound < lower_bound:
        return False
    else:
        return True

#Updates the top k alignments based on a new alignment
def updateTopAlignments(align):
    top_alignments.sort(key=lambda x: x.score)
    lowest_score = top_alignments[0].score
    if lowest_score < align.score:
        top_alignments[0] = align
    return continueDFS(align.correspVertices, align.score, lowest_score)


#Main GraphMatch recurrence (Fig 3)
def GraphMatch(W_in, W_prime_in, corr, G_prime, G0, G):
    global max_score
    for vi in corr:
        W = copy.copy(W_in)
        if vi not in W:
            W.append(vi)
            subgraph = G0.subgraph(W)
            if nx.is_connected(subgraph) :
                for vij in corr[vi]:
                    W_prime = copy.copy(W_prime_in) #Making copies because with recurrence, don't want to have pass by reference be an issue
                    W_prime.append(vij)
                    W_prime_str = []
                    for v in W_prime:
                        W_prime_str.append(v.stringifyVertex())
                    if isValidSolution(W, W_prime, G0, G_prime) and not T.hasPath(W_prime_str):
                        alignment = ScoreAlignment(W_prime, G, G0)
                        continue_search = updateTopAlignments(alignment)
                        if continue_search:
                            GraphMatch(W, W_prime, corr, G_prime, G0, G)
                    T.addPath(W_prime_str)

#Print the parameters selected by the user, and set up the specified parameters
def command_setup(args):
    print "Running GraphMatch...\n"
    print "Parameters chosen: "
    print "Query Graph: %s\nReference Graph: %s\nCorrespondences File: %s" % (args.q, args.i, args.c)
    print "Number of indels allowed: %d\nIndel Penalty: %.2f\nNon-associated query vertex penalty: %.2f\nNumber of top results: %d" % (args.m, args.ip, args.np, args.k)
    if args.s:
        print "Scoring Method: negative log(E-value)"
    else:
        print "Scoring Method: Raw match score"
        print "WARNING: due to scores being positive or negative for this approach, no branch and bound available."
    print
    #Setting up the parameters
    global m
    global indel_penalty
    global non_assoc_vertex_penalty
    global score_method
    global top_alignments
    m = args.m
    indel_penalty = args.ip
    non_assoc_vertex_penalty = args.np
    score_method = args.s
    #Instantiate the top alignments with low scoring Alignments (So they will get replaced)
    for i in range(0, args.k):
        a = Alignment(float("-inf"), [], [], [], [])
        top_alignments.append(a)

#Printing top alignment results
def print_results(G_refGraph):
    result_count = 1
    for align in sorted(top_alignments, key=lambda x: x.score, reverse=True):
        print "Result %d - score: %f" % (result_count, align.score)
        print
        print "Matches: (query \ corresponding vertex)"
        list_nodes = []
        for w in align.correspVertices:
            print "%s \t %s" % (w.queryVertex, w.name)
            list_nodes.append(w.name)
        print
        print "Edges:"
        for i in range(0, len(align.query_edges)):
            print "%s\t--\t%s\t\t" % (align.query_edges[i][0], align.query_edges[i][1]),
            for j in range(0, len(align.ref_edges[i]) - 1):
                print "%s\t--\t" % align.ref_edges[i][j],
            print align.ref_edges[i][len(align.ref_edges[i]) - 1]

        #Adding indel vertices to a list for getting the G subgraph ready
        for i in align.indelVertices:
            list_nodes.append(i)
        result = nx.subgraph(G_refGraph, list_nodes)
        for n in result.node:
            if n in align.indelVertices:
                result.node[n]['color'] = 'red'
            else:
                result.node[n]['color'] = 'blue'
        output_str = "result" + str(result_count) + ".dot"
        write_dot(result, output_str)
        result_count += 1
        print


def main():
    # Reading in the arguments from the command line
    parser = argparse.ArgumentParser(description='GraphMatch')
    parser.add_argument('q', type=str, help='file name containing query graph')
    parser.add_argument('i', type=str, help='file name containing reference graph')
    parser.add_argument('c', type=str, help='file name containing correspondences between vertices')
    parser.add_argument('-m', type=int, help = 'Number of indels allowed in resulting subgraph [default=1]', default=1)
    parser.add_argument('-ip', type=float, help='Indel Penalty [default=0.1]', default=0.1)
    parser.add_argument('-np', type=float, help='Penalty for query vertex being missing from resulting subgraph [default=1]', default=1)
    parser.add_argument('-s', type=int, help='Match score form: 0-raw match score; 1- -log(match score) [default=1]', default=1)
    parser.add_argument('-k', type=int, help='Number of top alignments to output [default=3]', default=3)
    args = parser.parse_args()

    command_setup(args)

    #Reading in graphs and correspondences
    G_refGraph = readGraph(args.i) #Node: String of corresponding vertex name
    G0_queryGraph = readGraph(args.q) #Node: String of query vertex name
    global num_query_nodes
    num_query_nodes = nx.number_of_nodes(G0_queryGraph)
    corr = readCorrespondances(args.c) #Dictionary: key=query vertex name; value=list of Corr_Vertex

    #Making G' graph in preparation for GraphMatch recurrence
    Gprime_graph = makeGraph_prime(corr, G0_queryGraph, G_refGraph) #Node: Stringified Corr_vertex

    GraphMatch([], [], corr, Gprime_graph, G0_queryGraph, G_refGraph) #Apply Graph Match with empty sets at first

    print_results(G_refGraph)


if __name__ == '__main__':
    main()