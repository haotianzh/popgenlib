from collections import OrderedDict, defaultdict
from .node import Node

class BaseTree(object):
    def __init__(self):
        self._root = None
        self._nodes = OrderedDict()
        self._hierarchical_tree = {}

    def __contains__(self, node):
        # if a node is in this tree
        if isinstance(node, Node):
            node = node.identifier
        if node in self._nodes:
            return True
        else:
            return False
    
    def __getitem__(self, identifier):
        # get a node from node list
        return self._nodes[identifier]

    def __len__(self):
        # obtain number of nodes in a tree
        return len(self._nodes)

    def from_newick(self, newick):
        # construct a tree from newick format string
        pass

    def create_node(self, name=None, identifier=None, parent=None):
        node = Node(name=name, identifier=identifier)
        self.add_node(node, parent)

    def add_node(self, node, parent=None):
        if node.identifier in self._nodes:
            raise Exception('cannot add the node that has already been in the tree.')
        if parent is None:
            if self._root:
                raise Exception('root has already existed and parent cannot be none.')
            self._root = node
            self._nodes[node.identifier] = node

            # /* set root and make its level as 0 */
            node.set_parent(None)
            return
        pid = parent.identifier if isinstance(parent, Node) else parent
        if not pid in self._nodes:
            raise Exception('parent not found in this tree.')
        self._nodes[node.identifier] = node

        # /* link node with its parent */
        node.set_parent(self[pid])
        self[pid].add_child(node)
        return

    def get_all_nodes(self, func=None):
        if func is None:
            return self._nodes
        return func(self._nodes)




