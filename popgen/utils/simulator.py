import msprime as msp
import numpy as np


class Simulate(object):
    """
    Using for creating simulated data.
    """
    def __init__(self):
        self.__simulator__ = '{name}/{version}'.format(name='msprime', version=msp.__version__)

    def __str__(self):
        return self.__simulator__


