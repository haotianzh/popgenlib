from time import time
import numpy as np
from popgen.utils.simulator import Simulator
from popgen.utils.treeutils import TraversalGenerator
from popgen.utils.statistics import cluster_ld, pairwise_ld
from popgen.utils.javautils import rentplus, rfdist
from popgen.utils.haputils import sliding_windows
start = time()

configs = {
    'sequence_length': 1e5,
    'population_size': 1e5,
    'rate': 1e-8
}
simulator = Simulator()
simulator.set_configs(configs)
ld = []
cl = []

for rep in simulator(100, 1):  # generate data using msprime
    print(rep.haplotype.nsites)

windows = sliding_windows(rep.haplotype, window_size=50)
res = rentplus(windows, num_thread=10)
rf = rfdist(res)





end = time()

print("running time %f" % (end - start))




