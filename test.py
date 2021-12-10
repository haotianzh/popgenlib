from time import time
from popgen.utils.simulator import Simulator
from popgen.utils.treeutils import TraversalGenerator
from popgen.utils.statistics import ld_matrix, ld_cluster
start = time()

configs = {
    'sequence_length': 1e5,
    'population_size': 1e5,
    'rate': 1e-9
}
simulator = Simulator()
simulator.set_configs(configs)
for sample in simulator(100, 1):
    ld = ld_matrix(sample.haplotype)
    cl = ld_cluster(ld)
    print(ld.shape, cl.shape)
import numpy as np


end = time()

print("running time %f" % (end - start))
