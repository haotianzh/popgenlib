from collections import OrderedDict
from uuid import uuid1

class Node(object):
    def __init__(self, name=None, identifier=None):
        self.name = name
        self.identifier = None
        if name is None:
            self.name = self.identifier

    @property
    def identifier(self):
        return self.identifier

    @identifier.setter
    def identifier(self, nid):
        _identifier = nid
        if nid is None:
            _identifier = str(uuid1())
        self.identifier = _identifier

