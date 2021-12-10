from popgen.base import BaseTree, Node
import pptree
from time import time
import numpy as np
from popgen.utils.treeutils import TraversalGenerator
from popgen.utils import from_newick
import msprime as msp
from popgen.utils.simulator import Simulator

start = time()

configs = {
    'sequence_length': 1e5,
    'population_size': 1e5,
    'rate': [1e-9, 1e-7]
}
simulator = Simulator()
simulator.set_configs(configs)
print(simulator.configs)
for sample in simulator(50, 5):
    print(sample.configs,'|', sample)


end = time()

print("running time %f" % (end - start))
