import numpy as np
from minisom import MiniSom  # Assuming you're using MiniSom for SOM implementation
import matplotlib.pyplot as plt
from Gene_Helper import GeneDataManager
import pandas as pd

# Assume the Gene_Helper class is properly defined and imported

class SOMClusterer:
    @staticmethod
    def perform_som_clustering(gene_manager, selected_features, grid_size=(10, 10), iterations=1000):
        gene_data_list = gene_manager.fetch_gene_data()

        data = []
        for gene_data in gene_data_list:
            feature_vector = [getattr(gene_data, feature) for feature in selected_features]
            data.append(feature_vector)

        data = np.array(data)
        data = (data - data.min(axis=0)) / (data.max(axis=0) - data.min(axis=0))

        som = MiniSom(grid_size[0], grid_size[1], data.shape[1], sigma=1.0, learning_rate=0.5)
        som.random_weights_init(data)
        som.train_batch(data, iterations, verbose=True)

        plt.figure(figsize=(10, 10))
        plt.pcolor(som.distance_map().T, cmap='bone_r')
        plt.colorbar()

        winning_neurons = np.array([som.winner(data[i]) for i in range(len(data))])

        cluster_data = {key: [] for key in [(i, j) for i in range(grid_size[0]) for j in range(grid_size[1])]}

        for i, j in cluster_data.keys():
            indices = np.where(np.logical_and(winning_neurons[:, 0] == i, winning_neurons[:, 1] == j))

            if len(indices[0]) > 0:
                cluster_data[(i, j)] = indices[0]

        for (i, j), indices in cluster_data.items():
            cluster_name = f'Cluster ({i}, {j})'
            print(f'\n{cluster_name}:')

            if len(indices) > 0:
                cluster_points = data[indices]
                feature_ranges = np.ptp(cluster_points, axis=0)

                for feature, feature_range in zip(selected_features, feature_ranges):
                    feature_min = np.min(cluster_points[:, selected_features.index(feature)])
                    feature_max = np.max(cluster_points[:, selected_features.index(feature)])
                    print(f'{feature} range: ({feature_min:.12f} - {feature_max:.12f})')

                print(f'Total Data Points in {cluster_name}: {len(indices)}')

        plt.show()

# Usage
#gene_manager = GeneDataManager('localhost', 'root', '720027924', 'pineal_data')
#selected_features = ['embryonic_day_21', 'postnatal_day_5', 'postnatal_day_20', 'postnatal_day_40',
                     #'embryonic_night_21', 'postnatal_night_5', 'postnatal_night_20', 'postnatal_night_40']

#SOMClusterer.perform_som_clustering(gene_manager, selected_features)

# Assuming gene_manager is an instance of GeneDataManager
gene_manager = GeneDataManager('localhost', 'root', '720027924', 'pineal_data')

# Get the significant expression data using the merged method
significant_expression = gene_manager.gene_data_manager_merged()

# Analyze the significant data
gene_manager.analyze_significant_data(significant_expression, gene_manager.fetch_gene_data())


