import networkx as nx
import pydotplus
from networkx.drawing.nx_pydot import write_dot

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

print nx.is_connected(G)
print G.has_edge(T, V)
G.add_edge(T, V)
print G.has_edge(T, V)

G.add_edge(V, U)
G.add_edge(V, Y)
G.add_edge(Y, H)

print nx.is_connected(G)

L = nx.shortest_path(G, T, U)
print "Length of shortest path: %d" % (len(L)-1) #Subtract 1 because both source and end node included in path

for l in L:
    print l

write_dot(G, "test.dot")

#Testing 'induced subgraphs'
X = G.subgraph([T, V, U])
write_dot(X, "test1.dot")
print X.has_edge(T, V)