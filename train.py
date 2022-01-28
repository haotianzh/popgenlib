from time import time
import numpy as np
import popgen
from popgen.utils.simulator import Simulator
from popgen.utils.treeutils import TraversalGenerator
from popgen.utils.statistics import cluster_ld, pairwise_ld, stat


from popgen.utils.javautils import rentplus, rfdist
from popgen.utils.utils import sliding_windows
# from popgen.nn.predict import predict_keras

# define simulation settings
configs = {
    'sequence_length': 1e5,
    'population_size': 1e5,
    'rate': 1e-8,
    'recombination_rate': 8.3e-8,
    'ploidy': 1
}

train_set = []
test_set = []

# create a simulator and cut them into small pieces
simulator = Simulator()
simulator.set_configs(configs)
for replicate in simulator(100, 10):
    simulator.update({'recombination_rate': np.random.rand() * 1e-7})
    rep = popgen.utils.filter_replicate(replicate)
    pieces = popgen.utils.cut_replicate(rep)
    train_set.extend([val.haplotype for val in pieces])

gens = popgen.utils.javautils.rentplus(train_set, num_thread=10)
print(gens)