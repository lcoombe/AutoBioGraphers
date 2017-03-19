#GraphMatch Code

import networkx as nx
import argparse
import re
from Corr_Vertex import Corr_Vertex
from Tree import Tree
from networkx.drawing.nx_pydot import write_dot
import copy
from math import log

#EXAMPLE USAGE: python GraphMatch.py test-data/query test-data/input test-data/corr_short

#GLOBAL VARIABLES

#Setting parameters based on Fig. 8 Caption from paper
non_assoc_vertex_penalty = 1.0
indel_penalty = 0.1
m = 1
#Tree to keep track of already enumerated subgraphs
T = Tree()
results = []
scores = []
score_method = 1

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

#     //Load in correspondances between vertices in G and G0
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

#     //Make G' --> will be used when enumerating solution graphs G0'
# V' = Union of all corresponding vertices (from dictionary structure)
# for each pair of corresponding vertices vij and vkl:
# 	If the vertices are not associated with the same vi, add edge (vij, vkl) to E' if (vi, vk) have edge in G0, and length of shortest path from vij to vkl in G is at most m+1
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

#TODO: Implement Scoring
#Scoring the alignment between two graphs G0 and G0'
#Potentially return the G0' subgraph? (including the indel vertices)
#Just alter parameters as needed -- I just put in placeholders
def ScoreAlignment(W, W_prime, G, G0):
# - For each vertex in G0', find the vertex correspondance in V0 (if exists)
# - For each vertex correspondance, add to the score
# - For each vertex in V0 that wasn't in a correspondance (as above), penalize by delta.
# - For each pair of associations (vi, vij), (vk, vkl) (such that vi and vj are in V0, and their corresponding vertices vij, vkl are in V0'):
# 	- If there is a (vi, vk) edge in G0:
# 		- Find the path length between vij and vkl in G0'
# 		- Penalize the score by -(path length-1)*delta

    score = 0
    assoc_vertices = []
    if score_method: #Use -log(corresp score) as match value. want to maximize this score.
        for node in G0.node:
            corresp = getCorrespondingVertex(node, W_prime)
            if corresp is None:
                score = score - non_assoc_vertex_penalty
            else:
                assoc_vertices.append(node)
                if score_method:
                    if corresp.score <= 1e-200:
                        score += -1.0*log(1e-200)
                    else:
                        score += -1.0*log(corresp.score)
                else:
                    score += corresp.score
        for i in range(0, len(assoc_vertices)):
            for j in range(i+1, len(assoc_vertices)):
                vi = assoc_vertices[i]
                vk = assoc_vertices[j]
                if G0.has_edge(vi, vk):
                    vij = getCorrespondingVertex(vi, W_prime)
                    vkl = getCorrespondingVertex(vk, W_prime)
                    path_G = nx.shortest_path(G, vij.name, vkl.name)
                    pathLen = len(path_G) - 2
                    score = score - pathLen*indel_penalty
    else: #Using raw correspondance score for matching -- want to maximize this score.
        for node in G0.node:
            corresp = getCorrespondingVertex(node, W_prime)
            if corresp is None:
                score = score - non_assoc_vertex_penalty
            else:
                assoc_vertices.append(node)
                score -= corresp.score
        for i in range(0, len(assoc_vertices)):
            for j in range(i+1, len(assoc_vertices)):
                vi = assoc_vertices[i]
                vk = assoc_vertices[j]
                if G0.has_edge(vi, vk):
                    vij = getCorrespondingVertex(vi, W_prime)
                    vkl = getCorrespondingVertex(vk, W_prime)
                    path_G = nx.shortest_path(G, vij.name, vkl.name)
                    pathLen = len(path_G) - 2
                    score = score - pathLen*indel_penalty
    return score



#Determining if the vertex combos are a valid solution
def isValidSolution(V0_plus, W_prime, G0, G_prime):
	# For each pair of  vertices (vi, vk) in V0+:
    for vi in V0_plus:
        for vk in V0_plus:
            if vi != vk and G0.has_edge(vi, vk):
                vij = getCorrespondingVertex(vi, W_prime)
                vkl = getCorrespondingVertex(vk, W_prime)
                if not G_prime.has_edge(vij.stringifyVertex(), vkl.stringifyVertex()):
                    return False
    return True


#Main GraphMatch recurrence (Fig 3)
def GraphMatch(W_in, W_prime_in, corr, G_prime, G0, G):
    for vi in corr:
        W = copy.copy(W_in)
        if vi not in W:
            W.append(vi) 	# V0+ = Union of W and vi (that set of vertices)
	# 		If the induced subgraph by adding vi to W is connected, and not already found (in T):
            subgraph = G0.subgraph(W)
            if nx.is_connected(subgraph) :
                for vij in corr[vi]:
                    W_prime = copy.copy(W_prime_in) #Making copies because with recurrence, don't want to have pass by reference be an issue
                    W_prime.append(vij)
                    W_prime_str = []
                    for v in W_prime:
                        W_prime_str.append(v.stringifyVertex())
                    if isValidSolution(W, W_prime, G0, G_prime) and not T.hasPath(W_prime_str):
                        score = ScoreAlignment(W, W_prime, G, G0) #TODO: Record the alignment and its score
                        results.append((W, W_prime)) #TODO: Replace this with keeping track of top k alignments
                        for w in W:
                            p = getCorrespondingVertex(w, W_prime)
                            print w + "\t" + p.name
                        print score
                        print
                        scores.append(score)
                        GraphMatch(W, W_prime, corr, G_prime, G0, G)
                    T.addPath(W_prime_str)

#NOTES:
#I'm having trouble with keeping track of the V0+ vertices in the tree.
#I think I'm seeing that if I add W to the path of T, then a lot of the combinations wanted aren't being generated
#Ie. If have vi = H, and adding S. Say both H and s have corresponding vertices: H:[a, b] and S:[c, d]
# Then get GraphMatch([H,S], [a,c], ..) call
# Eventually, that will return and add [H, S] to the Tree.
# Then, GraphMatch([H,S], [a,d], ..) might be called, but then it would appear to have that path in T already, and
# wouldn't go forward even though this could be a valid solution

#For now, just adding the W_prime paths, since we know we don't want to look at multiple of those collections of vertices multiple times.

def main():
    # Reading in the arguments from the command line
    parser = argparse.ArgumentParser(description='GraphMatch')
    parser.add_argument('q', type=str, help='file name containing query graph')
    parser.add_argument('i', type=str, help='file name containing reference graph')
    parser.add_argument('c', type=str, help='file name containing correspondences between vertices')
    parser.add_argument('-m', type=int, help = 'Number of indels allowed in resulting subgraph [default=1]', default=1)
    parser.add_argument('-ip', type=float, help='Indel Penalty [default=0.1]', default=-0.1)
    parser.add_argument('-np', type=float, help='Penalty for query vertex being missing from resulting subgraph [default=1]', default=-1)
    parser.add_argument('-s', type=int, help='Match score form: 0-raw match score; 1- -log(match score)', default=1)
    args = parser.parse_args()

    global m
    global indel_penalty
    global non_assoc_vertex_penalty
    global score_method
    m = args.m
    indel_penalty = args.ip
    non_assoc_vertex_penalty = args.np
    score_method = args.s


    G_refGraph = readGraph(args.i) #Node: String of corresponding vertex name
    G0_queryGraph = readGraph(args.q) #Node: String of query vertex name

    corr = readCorrespondances(args.c) #Dictionary: key=query vertex name; value=list of Corr_Vertex
    Gprime_graph = makeGraph_prime(corr, G0_queryGraph, G_refGraph) #Node: Stringified Corr_vertex (See class function)

    GraphMatch([], [], corr, Gprime_graph, G0_queryGraph, G_refGraph) #Apply Graph Match with empty sets at first

    write_dot(G0_queryGraph, 'G0.dot')
    write_dot(Gprime_graph, 'G_prime.dot')

#TEMPORARY: PRINTING RESULTS (TODO: Keep track of top k alignments)
    # for r in results:
    #     W = r[0]
    #     W_prime = r[1]
    #     for w in W:
    #         p = getCorrespondingVertex(w, W_prime)
    #         print w + "\t" + p.name
    #     print

    print sorted(scores)

if __name__ == '__main__':
    main()