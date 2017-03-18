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
results = []

T1 = Tree()
ll = []

prefix_list = []

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


#TODO: Implement Scoring
#Scoring the alignment between two graphs G0 and G0'
#Potentially return the G0' subgraph? (including the indel vertices)
#Just alter parameters as needed -- I just put in placeholders
def ScoreAlignment(W, W_prime, G, G0):
    #print "TODO: Score Alignment code"
    x = 1
# - For each vertex in G0', find the vertex correspondance in V0 (if exists)
# - For each vertex correspondance, add to the score
# - For each vertex in V0 that wasn't in a correspondance (as above), penalize by delta.
# - For each pair of associations (vi, vij), (vk, vkl) (such that vi and vj are in V0, and their corresponding vertices vij, vkl are in V0'):
# 	- If there is a (vi, vk) edge in G0:
# 		- Find the path length between vij and vkl in G0'
# 		- Penalize the score by -(path length-1)*delta

#Given a vertex vi, and a list of corresponding vertex, return the corresponding vertex of vi
def getCorrespondingVertex(v, correspList):
    for c in correspList:
        if c.queryVertex == v:
            return c
    print "FAIL! in getCorrespondingVertex" #Shouldn't get here, as should have a corresponding vertex present

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

def areListsSame(list1, list2):
    if len(list1) != len(list2):
        return False
    for i in range(0, len(list1)):
        if list1[i] != list2[i]:
            return False
    return True

def checkLists(list1):

    for entry in prefix_list:
        if areListsSame(list1, entry):
            return True

    return False

#Main GraphMatch recurrence (Fig 3)
def GraphMatch(W_in, W_prime_in, corr, G_prime, G0):
    for vi in corr:
        W = copy.copy(W_in)
        if vi not in W:
            W.append(vi) 	# V0+ = Union of W and vi (that set of vertices)
	# 		If the induced subgraph by adding vi to W is connected, and not already found (in T):
            subgraph = G0.subgraph(W)
            if T.hasPath(W):
                same_prefix = checkLists(W) #Saving those paths with the same sequence of W --> meaning that they are still valid to check.

            if nx.is_connected(subgraph) and ((not T.hasPath(W)) or same_prefix):
                count = 0
                for vij in corr[vi]:
                    count += 1
                    W_prime = copy.copy(W_prime_in) #Making copies because with recurrence, don't want to have pass by reference be an issue
                    W_prime.append(vij)
                    if isValidSolution(W, W_prime, G0, G_prime):
                        #score = ScoreAlignment(W, W_prime, G, G0) #TODO: Record the alignment and its score

                        results.append((W, W_prime)) #TODO: Replace this with keeping track of top k alignments
                        GraphMatch(W, W_prime, corr, G_prime, G0)
                T.addPath(W)
                prefix_list.append(W) #Need this to rescue the valid solutions that are otherwise lost when calling GraphMatch recursion multiple times with same W, different W'

#NOTES:
#Bit hacky, but one way I figured we could rescue the solutions is to keep track of a list of what we added to the tree
# What we care about is if there is exactly the same list, with W in the exact same order in prefix_list
# Then, can assume this is a GraphMatch run with same W, but different W_prime
# It's a bit weird - open to suggestions for better ways to do it

def main():
    # Reading in the arguments from the command line
    parser = argparse.ArgumentParser(description='GraphMatch')
    parser.add_argument('q', type=str, help='file name containing query graph')
    parser.add_argument('i', type=str, help='file name containing reference graph')
    parser.add_argument('c', type=str, help='file name containing correspondences between vertices')
    args = parser.parse_args()

    G_refGraph = readGraph(args.i) #Node: String of corresponding vertex name
    G0_queryGraph = readGraph(args.q) #Node: String of query vertex name

    corr = readCorrespondances(args.c) #Dictionary: key=query vertex name; value=list of Corr_Vertex
    Gprime_graph = makeGraph_prime(corr, G0_queryGraph, G_refGraph) #Node: Stringified Corr_vertex (See class function)

    GraphMatch([], [], corr, Gprime_graph, G0_queryGraph) #Apply Graph Match with empty sets at first

    write_dot(G0_queryGraph, 'G0.dot')
    write_dot(Gprime_graph, 'G_prime.dot')

#TEMPORARY: PRINTING RESULTS (TODO: Keep track of top k alignments)
    for r in results:
        W = r[0]
        W_prime = r[1]
        for w in sorted(W):
            p = getCorrespondingVertex(w, W_prime)
            print w + "\t" + p.name
        print



if __name__ == '__main__':
    main()