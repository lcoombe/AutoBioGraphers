import networkx as nx
import pydotplus
from networkx.drawing.nx_pydot import write_dot

#NOTE: All this testing is based on playing around with the Networkx documentation, which can be found at https://networkx.readthedocs.io/en/stable/

class Testing: #Because for GraphMatch, we might represent nodes using a class

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return "%d, %d" % (self.x, self.y)

    def add(self):
        return self.x + self.y


G = nx.Graph()

T = Testing(1, 2)
V = Testing(3, 4)
Y = Testing(3, 5)
H = Testing(6, 7)
U = Testing(2, 1)

G.add_node(T)
G.add_node(V)
G.add_node(Y)
G.add_node(H)
G.add_node(U)
G.add_node(U)

print "Printing out the nodes of graph G:"
for node in G.node:
	print node 
print 


print "Is G connected?: %s" % nx.is_connected(G)
print "Does G have an edge from nodes T and V? %s" % G.has_edge(T, V)
print "Does G have an edge from nodes V to T? %s" % G.has_edge(V, T)
print "Now adding edge from T and V"
G.add_edge(T, V)
print "Now, does G have an edge from nodes T and V? %s" % G.has_edge(T, V)
print "Now, does G have an edge from nodes V to T? %s" % G.has_edge(V, T)

print "Adding more edges..."
G.add_edge(V, U)
G.add_edge(V, Y)
G.add_edge(Y, H)
G.add_edge(Testing(3, 4), Testing(10, 19))

print "Showing all nodes: "
for l in G.node:
	print l 

print "Is G connected? %s" % nx.is_connected(G)

L = nx.shortest_path(G, T, U)
print "Length of shortest path from T to U: %d" % (len(L)-1) #Subtract 1 because both source and end node included in path

print "The nodes in the shortest path:"
for l in L:
    print l

#Writing out G to dot file
write_dot(G, "test.dot")

#Testing 'induced subgraphs'
X = G.subgraph([T, V, U])
write_dot(X, "test1.dot")
print "Does X (induced subgraph of G with nodes T, V, U) have an edge from T to V? %s" % X.has_edge(T, V)