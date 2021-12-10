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
    'population_size': 1e5
}
simulator = Simulator()
simulator.set_configs(configs)
for sample in simulator(50, 5):
    simulator.update({'rate': np.random.rand() * 1e-9})
    print(sample.configs['rate'],'|', sample)


end = time()

print("running time %f" % (end - start))
