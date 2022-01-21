from time import time
import numpy as np
import popgen
from popgen.utils.simulator import Simulator
from popgen.utils.treeutils import TraversalGenerator
from popgen.utils.statistics import cluster_ld, pairwise_ld, stat
from popgen.utils.javautils import rentplus, rfdist
from popgen.utils.utils import sliding_windows
from popgen.nn.predict import predict_keras

start = time()

configs = {
    'sequence_length': 1e5,
    'population_size': 1e5,
    'rate': 1e-8,
    'recombination_rate': 3.4e-8,
    'ploidy': 1
}

simulator = Simulator()
simulator.set_configs(configs)

for rep in simulator(100, 1):  # generate data using msprime
    print(rep.haplotype)
print(rep.haplotype.nsites)
# cut haplotype
filtered_hap = popgen.utils.utils.filter_none_mutation(rep.haplotype)
windows2 = sliding_windows(filtered_hap, window_size=50)
# cut replicate
filtered_rep = popgen.utils.utils.filter_replicate(rep)
print(filtered_rep.ts)
replicates = popgen.utils.utils.cut_replicate(filtered_rep, window_size=50)
print(replicates[0].haplotype.matrix.shape)
windows1 = []
for win in replicates:
    windows1.append(win.haplotype)


def predict(windows):
    res = rentplus(windows, num_thread=20)
    rf = rfdist(res)
    lds = []
    cls = []
    for win in windows:
        ld = pairwise_ld(win)
        cl = cluster_ld(ld)
        lds.append(ld)
        cls.append(cl)
    rf = np.array(rf)
    lds = np.array(lds)
    cls = np.array(cls)
    features = np.concatenate([lds[..., np.newaxis], cls[..., np.newaxis], rf[..., np.newaxis]/2/97], axis=-1).astype(np.float16)
    rates = predict_keras(features, model_path='D:/Py/popgen/popgen/nn_models/snp50_rho200.mdl')
    scaled_rates, a, b = stat(rates, rep.haplotype.positions, window_size=50, step_size=50, bin_width=1e4, sequence_length=configs['sequence_length'], ne=configs['population_size'])
    print(scaled_rates.mean())


predict(windows1)
predict(windows2)


#
# print((rep.haplotype.matrix == 2).sum())
#
# print("running time %f" % (end - start))




