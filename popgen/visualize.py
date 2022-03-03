import popgen
import numpy as np
import matplotlib.pyplot as plt

n_sam = 50
n_rep = 100
window_size = 100

configs = {
    'sequence_length': 1e5,
    'population_size': 1e5,
    'rate': 1e-8,
    'recombination_rate': 8.3e-8,
    'ploidy': 1
}

recombination_rates = [1e-10, 1e-9, 1e-8, 5e-8, 1e-7, 5e-7]

def transform(newicks):
    trees = [popgen.utils.treeutils.from_newick(newick) for newick in newicks]
    genealogies = []
    union = set()
    for tree in trees:
        splits = tree.get_splits()
        genealogies.append(splits)
        union.update(splits)
    union_list = list(union)
    union_list = sorted(union_list, key=lambda x: len(x), reverse=True)
    map_tensor = np.zeros([len(union), window_size])
    for idx, split in enumerate(union_list):
        for idy, genealogy in enumerate(genealogies):
            if split in genealogy:
                map_tensor[idx, idy] = 1
    return map_tensor

data = []
haplotypes = []
rates = []
simulator = popgen.Simulator()
simulator.set_configs(configs)
for i, rate in enumerate(recombination_rates):
    simulator.update({'recombination_rate': rate})
    for replicate in simulator(n_sam, 1):
        rep = popgen.utils.filter_replicate(replicate)
        # print(rep.ts.rr, rep.ts.mr, rep.haplotype.positions[-1]-rep.haplotype.positions[0])
        pieces = popgen.utils.cut_replicate(rep, window_size=window_size)
        data.extend([(val.haplotype, val.configs['recombination_rate']) for val in pieces])

for d in data:
    haplotype, rate = d
    haplotypes.append(haplotype)
    rates.append(rate)

gens = popgen.utils.rentplus(haplotypes, num_thread=10)
gens = [[str(val) for val in z] for z in gens]

def plot(newicks, rate):
    map_tensor = transform(newicks)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(rate)
    cax = ax.matshow(map_tensor)
    fig.colorbar(cax)
    plt.show()