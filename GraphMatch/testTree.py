__author__ = 'laurencoombe'

from Tree import Tree

T = Tree()
T.addPath(["a", 'b', 'bac'])
T.addPath(['t', 'e', 'f'])
T.addPath(['a', 'b'])
T.addPath(['a', 'b'])
T.addPath(['b', 'c'])
T.addPath(['a', 'b', 'x', 'y'])

print T.hasPath(['b', 'x', 'a', 'y']) #Should be true
print T.hasPath(['e', 'f']) #Should be false, as no path was added consisting of only these nodes

T.printTree()
T.addPath(['f', 'e'])

print T.hasPath(['e', 'f']) #Now should be True, as we added this exact path

print
T.printTree()