from time import time
import numpy as np
from sklearn.metrics import mean_absolute_error



import popgen
from popgen.utils.simulator import Simulator
from popgen.utils.treeutils import TraversalGenerator
from popgen.utils.statistics import cluster_ld, pairwise_ld, stat
from popgen.utils.javautils import rentplus, rfdist
from popgen.utils.utils import sliding_windows




import torch
from torch import nn
from torch.utils.data import DataLoader
from popgen.nn.network import SplitNet

# define simulation settings
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


class GenDataSet(torch.utils.data.Dataset):
    def __init__(self, train_set, transform=None):
        self.transform = transform
        self.haplotypes = []
        self.rates = []
        for hap, rate in train_set:
            self.haplotypes.append(hap)
            self.rates.append(rate)
        genealogies = popgen.utils.javautils.rentplus(self.haplotypes, num_thread=10)
        self.genealogies = [[str(v) for v in val] for val in genealogies]

    def __len__(self):
        return len(self.haplotypes)

    def __getitem__(self, idx):
        newicks = self.genealogies[idx]
        rate = self.rates[idx]
        # print(newicks)
        distance = self.haplotypes[idx].positions[-1] - self.haplotypes[idx].positions[0]
        t1, t2 = self.transform(newicks)
        return t1, t2, rate, distance


def transform(newicks):
    trees = [popgen.utils.treeutils.from_newick(newick) for newick in newicks]
    genealogies = []
    union = set()
    for tree in trees:
        splits = tree.get_splits()
        genealogies.append(splits)
        union.update(splits)
    splits_tensor = torch.zeros([len(union), n_sam])
    map_tensor = torch.zeros([len(union), window_size])
    # print(splits_tensor.shape, map_tensor.shape)
    for id_split, split in enumerate(union):
        for leaf in split:
            splits_tensor[id_split, int(leaf) - 1] = 1.0
        for j, genealogy in enumerate(genealogies):
            if split in genealogy:
                map_tensor[id_split][j] = 1.0
    return splits_tensor, map_tensor


def train():
    training_history = []

    def validate(model, testdata):
        predicted_rates = []
        truth = []
        for d in testdata:
            splits, maps, rate, distance = d
            splits, maps, rate = splits.to(device), maps.to(device), torch.tensor(rate).to(device)
            with torch.no_grad():
                output = model(splits, maps)
            psr = output.item() / configs['population_size'] / 2 / configs['ploidy'] / distance
            predicted_rates.append(psr)
            truth.append(rate)
        mae = mean_absolute_error(truth, predicted_rates)
        training_history.append(mae)
        print(f'truth: {mae}')

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    data = []
    # create a simulator and cut them into small pieces
    # !!!!!!!!!!!!!!!  have issues here
    simulator = Simulator()
    simulator.set_configs(configs)
    for i in range(n_rep):
        simulator.update({'recombination_rate': np.random.rand() * 1e-7})
        for replicate in simulator(n_sam, 1):
            rep = popgen.utils.filter_replicate(replicate)
            print(rep.ts.rr, rep.ts.mr)
            pieces = popgen.utils.cut_replicate(rep, window_size=window_size)
            data.extend([(val.haplotype, val.configs['recombination_rate']) for val in pieces])
    train_ratio = 0.8
    train_set = data[:int(len(data)*train_ratio)]
    test_set = data[int(len(data)*train_ratio):]
    print(f'train_set_length: {len(train_set)}, test_set_length: {len(test_set)}')
    dataset = GenDataSet(train_set, transform=transform)
    dataset_test = GenDataSet(test_set, transform=transform)
    model = SplitNet(num_nodes=n_sam)
    print(f'training model: {model}')
    model.to(device=device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    # train the model
    epochs = 20
    loss_fn = torch.nn.MSELoss()
    for e in range(epochs):
        for i, d in enumerate(dataset):

            splits, maps, rate, distance = d
            splits, maps, rate = splits.to(device), maps.to(device), torch.tensor(rate).to(device)
            # print('split', splits.shape, 'maps', maps.shape)
            psr = rate * configs['population_size'] * 2 * configs['ploidy'] * distance
            # print('psr', psr)
            optimizer.zero_grad()
            output = model(splits, maps)
            loss = torch.abs(output - psr)
            loss.backward()
            optimizer.step()
            print(f'training loss: {loss.item()}')
        validate(model, dataset_test)
    return training_history

if __name__ == '__main__':
    history = train()
    print(history)