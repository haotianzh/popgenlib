from time import time
import numpy as np
import popgen
from popgen.utils.simulator import Simulator
from popgen.utils.treeutils import TraversalGenerator
from popgen.utils.statistics import cluster_ld, pairwise_ld, stat
from popgen.utils.javautils import rentplus, rfdist
from popgen.utils.utils import sliding_windows

# from popgen.nn.predict import predict_keras

start = time()

# configs = {
#     'sequence_length': 1e5,
#     'population_size': 1e5,
#     'rate': 1e-8,
#     'recombination_rate': 8.3e-8,
#     'ploidy': 1
# }
#
# simulator = Simulator()
# simulator.set_configs(configs)
#
# for rep in simulator(100, 1):  # generate data using msprime
#     print(rep.haplotype)
# print(rep.haplotype.nsites)
# # cut haplotype
# filtered_hap = popgen.utils.utils.filter_none_mutation(rep.haplotype)
# windows2 = sliding_windows(filtered_hap, window_size=50)
# # cut replicate
# filtered_rep = popgen.utils.utils.filter_replicate(rep)
# print(filtered_rep.ts)
# replicates = popgen.utils.utils.cut_replicate(filtered_rep, window_size=50)
# print(replicates[0].haplotype.matrix.shape)
# windows1 = []
# for win in replicates:
#     windows1.append(win.haplotype)
#
#
# def predict(windows):
#     res = rentplus(windows, num_thread=20)
#     rf = rfdist(res)
#     lds = []
#     cls = []
#     for win in windows:
#         ld = pairwise_ld(win)
#         cl = cluster_ld(ld)
#         lds.append(ld)
#         cls.append(cl)
#     rf = np.array(rf)
#     lds = np.array(lds)
#     cls = np.array(cls)
#     features = np.concatenate([lds[..., np.newaxis], cls[..., np.newaxis], rf[..., np.newaxis]/2/97], axis=-1).astype(np.float16)
#     rates = predict_keras(features, model_path='D:/Py/popgen/popgen/nn_models/snp50_rho200.mdl')
#     scaled_rates, a, b = stat(rates, rep.haplotype.positions, window_size=50, step_size=50, bin_width=1e4, sequence_length=configs['sequence_length'], ne=configs['population_size'])
#     print(scaled_rates.mean())
#
#
# predict(windows1)
# predict(windows2)
#
import popgen
import pandas as pd


def do_once(tree):
    tree = popgen.utils.treeutils.from_newick(tree)
    splits = tree.get_splits()
    sorted_splits = sorted(splits, key=lambda x: len(x))
    data = {}
    for split in sorted_splits:
        data[split] = []
        for leaf in tree.get_leaves():
            data[split].append(int(leaf in split))
    df = pd.DataFrame(data, index=tree.get_leaves())
    return sorted_splits

trees = '''(((((((((17,8),20),5),((11,14),6)),(((4,9),13),7)),(((10,12),2),(16,19))),1),(18,3)),15)
((((((((((17,8),20),5),((11,14),6)),(((4,9),13),7)),((16,19),1)),((10,12),2)),3),18),15)
((((((((((17,8),20),5),((11,14),6)),(((4,9),13),7)),(18,3)),2),((10,12),(16,19))),1),15)
(((((((((((17,8),20),5),(((4,9),13),7)),(((11,14),6),3)),18),2),(16,19)),(10,12)),1),15)
((((((((((17,8),20),5),(((4,9),13),7)),(((11,14),6),3)),18),2),((10,12),(16,19))),1),15)
((((((((((17,8),20),5),(((4,9),7),13)),18),(((11,14),6),3)),2),((10,12),(16,19))),1),15)
((((((((((10,12),(19,2)),((11,14),6)),3),1),(((17,8),20),5)),(16,18)),((4,9),13)),7),15)
(((((((((((12,19),10),2),((11,14),6)),3),1),(((17,8),20),5)),(16,18)),((4,9),13)),7),15)
((((((((((19,2),6),3),(11,14)),((10,12),1)),(((17,8),20),5)),(((4,9),13),18)),16),7),15)
(((((((((((17,8),20),5),3),13),((19,2),6)),((11,14),12)),18),16),15),((((10,4),9),7),1))
((((((((((17,8),20),5),(13,3)),(19,2)),((11,14),(12,6))),18),16),((((10,4),9),7),1)),15)
(((((((((17,8),20),5),(13,3)),(19,2)),11),(((((12,14),6),16),15),1)),7),((10,4),(18,9)))
(((((((((19,2),6),15),7),5),(((12,14),11),1)),(((17,8),20),(13,3))),16),((10,4),(18,9)))
(((((((15,6),(19,2)),7),5),((((17,8),20),13),3)),((((11,14),12),1),16)),((10,4),(18,9)))
(((((((((19,2),6),15),7),5),((((17,8),20),13),(11,3))),(((12,14),1),16)),(18,9)),(10,4))'''
trees = [tree+';' for tree in trees.split('\n')]
print(trees)
all_splits = []

for tree in trees:
    split = do_once(tree)
    all_splits.append(split)

#
# print((rep.haplotype.matrix == 2).sum())
#
# print("running time %f" % (end - start))
