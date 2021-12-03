from collections import OrderedDict
from uuid import uuid1
from copy import deepcopy


class Node(object):
    def __init__(self, name=None, identifier=None):
        self.identifier = identifier
        if name is None:
            self.name = self.identifier
        else:
            self.name = name
        self.parent = {}
        self.children = OrderedDict()
        self.childs = []


    def add_child(self, node):
        if node in self.children:
            raise Exception('node has already been added.')
        if node in self.parent:
            raise Exception('parent cannot be added as a child.')
        self.children[node.identifier] = node
        self.childs.append(node)

    def set_parent(self, parent):
        if parent is None:
            self.parent = {}
        else:
            self.parent[parent.identifier] = parent

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

if __name__ == '__main__':
    node = Node()
    print(node.identifier, node.name)