import popgen
import numpy as np
from popgen.utils.statistics import pairwise_ld, cluster_ld

configs = {
    'sequence_length': 5e5,
    'population_size': 1e5,
    'rate': 1e-8,
    'recombination_rate': 3.9e-8,
    'ploidy': 1
}
n_sams = 5 # for test only
n_pop = 100
n_per_draw = 5
global_window_size = 1000
window_size = 50
r_max = 1e-7
r_min = 1e-9
# return the recombination rate given index i
def recomb_rate(i, rmin=r_min, rmax=r_max):
    rate = np.exp((np.log(rmax)-np.log(rmin))/(n_sams-1)*i + np.log(rmin))
    return rate

haplotypes = []
genealogies = []
rhos = []
# init simulator
simulator = popgen.Simulator(configs)

# generate random haplotypes and infer their genealogies.
for i in range(n_sams):
    rr = recomb_rate(i)
    print(rr)
    simulator.update({'recombination_rate': rr})
    for data in simulator(n_pop, 1):
        data = popgen.utils.utils.filter_replicate(data)
        reps = popgen.utils.cut_replicate(data, window_size=global_window_size)
        haps = [rep.haplotype for rep in reps]
        inferred_genealogies = popgen.utils.rentplus(haps, num_thread=10)
        for hap, gen in zip(haps, inferred_genealogies):
            for j in range(n_per_draw):
                start = np.random.choice(range(hap.nsites - window_size))
                haplotypes.append(hap.matrix[:, start:start+window_size])
                genealogies.append(gen[start: start+window_size])
                length = hap.positions[start+window_size] - hap.positions[start]
                scaled_rho = 2*configs['population_size']*rr*length
                rhos.append(scaled_rho)

# build the whole set, train set and test set. compute ld, ld_cluster and rf_dist.

rfdistance = popgen.utils.rfdist([list(val) for val in genealogies])
lds = []
cls = []
for hap, gen, rf in zip(haplotypes, genealogies, rfdistance):
    ld = pairwise_ld(hap)
    cl = cluster_ld(ld)
    lds.append(ld)
    cls.append(cl)
lds = np.array(lds)[..., np.newaxis]
cls = np.array(cls)[..., np.newaxis]
rfs = np.array(rfdistance)[..., np.newaxis]
data = np.concatenate([lds,cls,rfs], axis=-1).astype(np.float16)


    
