from . import Node, BaseTree, Haplotype


class Replicate(object):

    def __init__(self, ts, configs):
        self.configs = configs
        self.ts = ts
        self.haplotype = Haplotype(ts)

    def genealogies(self, branch_length=False):
        genealogies = []
        for tree in self.ts.trees():
            genealogies.append(tree.newick(include_branch_lengths=branch_length))
        return genealogies

    def __str__(self):
        return "%d samples, %d sites, %d topologies" % \
               (self.haplotype.nsamples, self.haplotype.nsites, len(self.genealogies()))
