from ..base import Haplotype, Sample
import msprime as msp
import numpy as np
import warnings


class Simulator(object):
    """
    Using for creating simulated data.
    calling method:
        nsam: number of samples
        nreps: number of replicates
    configs:
        rate: mutation rate (per base)
        recombination_rate: recombination rate (per base)
        sequence_length: length of simulated genome (bp)
        population_size(or demography): either one of them should be specified
        ploidy: haploid or diploid
    """

    def __init__(self):
        self.__simulator__ = '{name}/{version}'.format(name='msprime', version=msp.__version__)
        self.configs = None
        self._mutation_configs = None
        self._ancestry_configs = None

    def __call__(self, nsam, nreps):
        """ Run simulation for nsam samples and nreps replicates """
        configs = self.configs
        if configs is None:
            raise Exception("no configuration for simulation found.")
        replicates = msp.sim_ancestry(samples=nsam, num_replicates=nreps, **self._ancestry_configs)
        # simulate mutations once for each of genealogical trees.
        for ts in replicates:
            mts = msp.sim_mutations(ts, **self._mutation_configs)
            sample = Sample(mts, self.configs)
            yield sample

    def update(self, u):
        assert isinstance(u, dict), Exception('key-value pairs must be provided.')
        configs = self.configs.copy()
        configs.update(u)
        self.set_configs(configs)

    def set_configs(self, configs):
        requires = ['sequence_length']
        if configs is None:
            raise Exception("configuration cannot be empty.")
        for key in requires:
            if key not in configs:
                raise Exception("%s must be set." % key)
        if 'population_size' in configs and 'demography' in configs:
            warnings.warn("both population size and demography detected, population size will be removed.")
            del (configs['population_size'])
        mutation_rate = configs['rate'] if 'rate' in configs else 1e-8
        self._mutation_configs = {'rate': mutation_rate}
        self._ancestry_configs = {'recombination_rate': 1e-8, 'ploidy': 1}
        self._ancestry_configs.update(configs)
        if self.configs is None:
            self.configs = dict()
        self.configs.update(self._mutation_configs)
        self.configs.update(self._ancestry_configs)
        if 'rate' in self._ancestry_configs:
            del (self._ancestry_configs['rate'])

    def __str__(self):
        return self.__simulator__

