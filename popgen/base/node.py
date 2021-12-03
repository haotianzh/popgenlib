from collections import OrderedDict
from uuid import uuid1
from copy import deepcopy


class Node(object):
    def __init__(self, name=None, identifier=None, branch=0):
        self.identifier = identifier
        if name is None:
            self.name = self.identifier
        else:
            self.name = name
        self.parent = {}
        self.branch = branch
        self.children = OrderedDict()
        self.children_list = []

    def is_root(self):
        return self.children


    def add_child(self, node):
        if node in self.children:
            raise Exception('node has already been added.')
        if node in self.parent:
            raise Exception('parent cannot be added as a child.')
        self.children[node.identifier] = node
        self.children_list.append(node)

    def set_parent(self, parent):
        if parent is None:
            self.parent = {}
            self.level = 0
        else:
            self.parent[parent.identifier] = parent
            self.level = parent.level + 1

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, nid):
        if nid is None:
            nid = str(uuid1())
        self._identifier = nid

    @property
    def level(self):
        return self._level

    @level.setter
    def level(self, value):
        if not isinstance(value, int):
            raise Exception('value of level should be an integer.')
        if value < 0:
            raise Exception('value cannot be negative.')
        self._level = value

    @property
    def branch(self):
        return self._branch
    
    @branch.setter
    def branch(self, value):
        if self.parent is None:
            raise Exception('node has no parent, branch cannot be assigned.')
        if isinstance(value, float) or isinstance(value, int):
            self._branch = value
        else:
            raise Exception('branch must be a number.')

    def set_branch(self, value):
        self.branch = value
    
if __name__ == '__main__':
    node = Node()
    print(node.identifier, node.name)