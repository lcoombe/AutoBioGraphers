#GraphMatch Code

import networkx as nx
import argparse
import re
from Corr_Vertex import Corr_Vertex
from Tree import Tree
from networkx.drawing.nx_pydot import write_dot
import copy

#EXAMPLE USAGE: python GraphMatch.py test-data/query test-data/input test-data/corr_short

#GLOBAL VARIABLES

#Setting parameters based on Fig. 8 Caption from paper
non_assoc_vertex_penalty = -1.0
indel_penalty = -0.1
m = 1
#Tree to keep track of already enumerated subgraphs
T = Tree()

#HELPER FUNCTIONS

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
# - for each V0, there is a set of correspondances, and a score
# - Possible format: Dictionary, where key is the V0, and value is the correspondance (small class)
#    * Perhaps create small class representing the 'corresponding vertex'
#    * Fields in the class: name, score, vertex corresponds to
#v1 v2 corr
#where v1 is a vertex in the query graph, v2 is a vertex in the input graph
def readCorrespondances(corrFileName):
    corr = {}
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
                    G_prime.add_edge(vij, vkl)
    return G_prime


#TODO: REQUIRED HELPER FUNCTIONS:

#Scoring the alignment between two graphs G0 and G0'
def ScoreAlignment(G0, G0_prime):
    print "TODO: Score Alignment code"
# - For each vertex in G0', find the vertex correspondance in V0 (if exists)
# - For each vertex correspondance, add to the score
# - For each vertex in V0 that wasn't in a correspondance (as above), penalize by delta.
# - For each pair of associations (vi, vij), (vk, vkl) (such that vi and vj are in V0, and their corresponding vertices vij, vkl are in V0'):
# 	- If there is a (vi, vk) edge in G0:
# 		- Find the path length between vij and vkl in G0'
# 		- Penalize the score by -(path length-1)*delta

#Determining if the vertex combos are a valid solution
def isValidSolution(V0_plus, W_prime):
    print "TODO: isValidSolution code"
	# For each pair of  vertices (vi, vk) in V0+:
	# 	if (vi, vk) is an edge in E0:
	# 		vij = corresponding vertex of vi in W'
	# 		vkl = corresponding vertex of vk in W'
	# 		(vij, vkl) must be an edge in E'
    return False

def GraphMatch(W_in, W_prime_in, corr, G_prime, G0):
    print "GraphMatch DFS goes here"

	# With the loop, we are doing this:
	# Enumerate all connected induced subgraphs of G0
	# (Induced means that we pick the vertices, and a pair of vertices are connected in the subgraph if they have an edge in G0)
	# - Each enumeration is a way to get V0+ (The vertices in G0 that have an association in G0')
	# - To avoid enumerating graphs multiple times, represent connected induced subgraph as a path in a tree where the vertices are in sorted order in the tree T, last vertex marked
    #
	# //W: Set of vertices in V0+ currently
	# //W': Set of vertices corresponding to V0+ vertices

	# 	for each vertex vi in V0, but not already in W:
    for vi in corr:
        W = copy.deepcopy(W_in)
        if vi not in W:
            W.append(vi) 	# V0+ = Union of W and vi (that set of vertices)
	# 		If the induced subgraph by adding vi to W is connected, and not already found (in T):
            subgraph = G0.subgraph(W)
            if not T.hasPath(W) and nx.is_connected(subgraph):
                for vij in corr[vi]:
                    W_prime = copy.deepcopy(W_prime_in) #Making copies because with recurrence, don't want to have pass by reference be an issue
                    W_prime.append(vij)
                    if isValidSolution(W, W_prime):
                        score = ScoreAlignment(G0, G_prime)
                        #TODO: Record the alignment and its score
                        GraphMatch(W, W_prime, G_prime, G0)
                T.addPath(W)


	# 			for each vertex vij in the correspondance list of vi:
	# 				Run subroutine 'isValidSolution(V0+, W' union vij)'
	# 				if valid solution:
	# 					Use ScoreAlignment subroutine to find the alignment score
	# 					record the score (Perhaps keep top 5 or something?)
	# 					GraphMatch(W union {vi}, W' union {vij})
	# 			add the W union {Vi} subgraph to T


def main():
    # Reading in the arguments from the command line
    parser = argparse.ArgumentParser(description='GraphMatch')
    parser.add_argument('q', type=str, help='file name containing query graph')
    parser.add_argument('i', type=str, help='file name containing reference graph')
    parser.add_argument('c', type=str, help='file name containing correspondences between vertices')
    args = parser.parse_args()

    G_refGraph = readGraph(args.i)
    G0_queryGraph = readGraph(args.q)

    corr = readCorrespondances(args.c)
    for vi in corr:
        print vi

    Gprime_graph = makeGraph_prime(corr, G0_queryGraph, G_refGraph)


    GraphMatch([], [], corr, Gprime_graph, G0_queryGraph) #Apply Graph Match with empty sets at first

if __name__ == '__main__':
    main()