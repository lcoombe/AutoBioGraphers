#PURPOSE OF CLASS: representing all the already tested induced subgraphs of G0.
# In the paper, a "*" is used to indicate the end of a path. For our purposes, then, a subgraph has already been tested if the path is in the graph,
# and the last node in the path has a "*", indicating this path exactly had been added to the graph previously.
# When adding to a graph, the nodes are sorted first. Then, the path is added, without adding duplicate nodes.
# If a new path being added has all nodes already in the tree, then the node corresponding to the last in the path has a "*" appended


class Tree:

    def __init__(self):
        self.root = {}

    def addPath(self, list_nodes):
        sorted_nodes = sorted(list_nodes)
        curr_parent = self.root
        num_nodes = len(sorted_nodes)
        count = 1
        for name in sorted_nodes:
            new_name = str(name) + "*"
            if str(name) not in curr_parent.keys():
                if new_name not in curr_parent.keys():
                    if count == num_nodes:
                        curr_parent[str(new_name)] = {} #Last in the sequence, not already there with or without *, so add node*
                        curr_parent = curr_parent[str(new_name)]
                    else:
                        curr_parent[str(name)] = {} #Not last in sequence, not already there with or without *, so add node
                        curr_parent = curr_parent[str(name)]
                else: #Key is in the node* form, not in the node form, but still already there
                    curr_parent = curr_parent[str(new_name)]
            else:
                #Key is there, in the node form
                if count == num_nodes: #The node is there, in node form, but we are at the last in the path, so convert to node*
                    curr_parent[new_name] = curr_parent.pop(str(name))
                    curr_parent = curr_parent[new_name]
                else:
                    curr_parent = curr_parent[str(name)] #We are not at the last in the path, so just update current parent

            count += 1


    def hasPath(self, list_nodes):
        sorted_nodes = sorted(list_nodes)
        curr_parent = self.root
        has_path = True
        num_nodes = len(sorted_nodes)
        count = 1
        for name in sorted_nodes:
            new_name = str(name) + "*"
            if str(name) not in curr_parent and new_name not in curr_parent: #Node with or without * is not in the level of nodes
                has_path = False
                break
            else:
                if count == num_nodes: #We are at the end of the sorted_nodes, so to be valid, this match has to have a "*"
                    if new_name not in curr_parent:
                        has_path = False
                        break
                else: #We are not at the end, so node* and node are both acceptable, but advance to the node that is in the directory
                    if str(name) in curr_parent:
                        curr_parent = curr_parent[str(name)]
                    else:
                        curr_parent = curr_parent[new_name]
                count += 1
        return has_path

    def printString(self, string, node, level):
        if node != {}:
            for n in node:
                string =  "- "*level + n
                print string
                self.printString(string, node[n], level+1)


    def printTree(self):
        self.printString("", self.root, 0)
