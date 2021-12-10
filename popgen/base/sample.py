from . import Node, BaseTree, Haplotype


class Sample(object):

    def __init__(self, ts, configs):
        self.configs = configs
        self.ts = ts
        self.haplotype = Haplotype(ts)
        self.genealogies = []
        for tree in ts.trees():
            self.genealogies.append(tree.newick())

    def __str__(self):
        return "%d samples, %d sites, %d topologies" % \
               (self.haplotype.nsamples, self.haplotype.nsites, len(self.genealogies))
