from tensorflow.keras.models import load_model
import msprime as msp
import popgen
from popgen.utils.statistics import *
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error
from scipy.stats import pearsonr


def to_ndarray(x):
    return np.array(x, dtype=np.float16)[..., np.newaxis]


rate_map = msp.RateMap(
    position=[0, 10000, 50000, 60000, 70000, 100000],
    rate = [1.4e-9, 7.5e-8, 1.2e-8, 3.2e-8, 1e-7]
)

CONFIGS = {
    'sequence_length': 1e5,
    'population_size': 1e5,
    'rate': 1e-8,
    'recombination_rate': rate_map,
    'ploidy': 1
}

simulator = popgen.Simulator(CONFIGS)
for data in simulator(100, 1):
    pass

data = popgen.utils.filter_replicate(data)

def calculate_metrics(rate_map, y):
    true_per_base_rates = np.zeros(CONFIGS['sequence_length'])
    for i in range(len(rate_map.position)-1):
        true_per_base_rates[rate_map.position[i]: rate_map.position[i+1]] = rate_map.rate[i]
    est_per_base_rates = np.zeros(CONFIGS['sequence_length'])
    for i in range(len(y)):
        start = data.haplotype.positions[i*50]
        end = data.haplotype.positions[(i+1)*50]
        est_per_base_rates[start:end] = y[i]
    mae = mean_absolute_error(true_per_base_rates, y_pred)





cut_data = popgen.utils.cut_replicate(data, window_size=1000)
haps = [rep.haplotype for rep in cut_data]
gens = popgen.utils.rentplus(haps)
all_gens = []
for gen in gens:
    for i in range(0,1000,50):
        all_gens.append(gen[i:i+50])
all_gens = [list(val) for val in all_gens]
print(all_gens)
rfs = popgen.utils.rfdist(all_gens)
lds = []
cls = []
recut_data = popgen.utils.cut_replicate(data, window_size=50)
for win in recut_data:
    ld = pairwise_ld(win.haplotype)
    cl = cluster_ld(ld)
    lds.append(ld)
    cls.append(cl)

lds = to_ndarray(lds)
cls = to_ndarray(cls)
rfs = to_ndarray(rfs)
lds1 = lds[:len(rfs)]
cls1 = cls[:len(rfs)]
print(lds.shape, cls.shape, rfs.shape)
x = np.concatenate([lds1, cls1, rfs/2/(100-3)], axis=-1)
model = load_model('models/model_epoch_19.hdf5')
y = model.predict(x)
scaled_rates, a, b = stat(y, data.haplotype.positions,
                           window_size=50, step_size=50, bin_width=1e4,
                           sequence_length=CONFIGS['sequence_length'],
                           ne=CONFIGS['population_size'])
print(scaled_rates)


# -------------- local gen comparison --------------

haps_recut = [val.haplotype for val in recut_data]
local_gens = popgen.utils.rentplus(haps_recut)
rfs_local = popgen.utils.rfdist(local_gens)
rfs_local = to_ndarray(rfs_local)
x2 = np.concatenate([lds, cls, rfs_local/2/(100-3)], axis=-1)
model_local = load_model('models/snp50_rho200.mdl')
y2 = model_local.predict(x2)
scaled_rates2, a2, b2 = stat(y2, data.haplotype.positions,
                           window_size=50, step_size=50, bin_width=1e4,
                           sequence_length=CONFIGS['sequence_length'],
                           ne=CONFIGS['population_size'])


# -------------- plot figures -----------------------
plt.plot(b[0],b[1], drawstyle='steps-mid', color='b')
plt.plot(b2[0],b2[1], drawstyle='steps-mid', color='r')
plt.show()