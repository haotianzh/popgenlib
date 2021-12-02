import numpy as np

class Hap(object):
    def __init__(self, ancestral_state=0):
        self.haplotype_matrix = None
        self.positions = None
        self.ancestral_state = ancestral_state

    @property
    def num_sites(self):
        if not self.haplotype_matrix:
            return self.haplotype_matrix.shape[0]
        else:
            return 0

    @property
    def num_samples(self):
        if not self.haplotype_matrix:
            return self.haplotype_matrix.shape[1]
        else:
            return 0

    def converse(self):
        # converts ancestral state from 0 to 1.
        pass

    def from_tskit(self, ts):
        # converts from tskit format.
        pass
