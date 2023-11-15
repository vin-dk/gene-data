import numpy as np
from minisom import MiniSom 
import matplotlib.pyplot as plt
from Gene_Helper import GeneDataManager
import pandas as pd
from sklearn.preprocessing import StandardScaler

x = 5
y = 5

class SOMClusterer:
    @staticmethod
    def perform_som_clustering(gene_manager, selected_features, grid_size=(x, y), iterations=1000):
        gene_data_list = gene_manager.fetch_gene_data()

        data = []
        for gene_data in gene_data_list:
            feature_vector = [getattr(gene_data, feature) for feature in selected_features]
            data.append(feature_vector)

        data = np.array(data)

        # Use StandardScaler for feature scaling
        scaler = StandardScaler()
        data = scaler.fit_transform(data)

        som = MiniSom(grid_size[0], grid_size[1], data.shape[1], sigma=1.0, learning_rate=0.5)
        som.random_weights_init(data)
        som.train_batch(data, iterations, verbose=True)

        winning_neurons = np.array([som.winner(data[i]) for i in range(len(data))])

        # Extract the actual grid size used by MiniSom
        actual_grid_size = som.get_weights().shape[:2]

        # Calculate cluster assignments for each data point based on winning neurons
        cluster_assignments = [i * actual_grid_size[1] + j for i, j in winning_neurons]

        # Display clusters
        for cluster_idx in range(actual_grid_size[0] * actual_grid_size[1]):
            cluster_name = f'Cluster {cluster_idx}'
            indices = np.where(np.array(cluster_assignments) == cluster_idx)[0]

            print(f'\n{cluster_name}:')

            if len(indices) > 0:
                cluster_points = data[indices]
                feature_ranges = np.ptp(cluster_points, axis=0)

                for feature, feature_range in zip(selected_features, feature_ranges):
                    feature_min = np.min(cluster_points[:, selected_features.index(feature)])
                    feature_max = np.max(cluster_points[:, selected_features.index(feature)])
                    print(f'{feature} range: ({feature_min:.12f} - {feature_max:.12f})')

                print(f'Total Data Points in {cluster_name}: {len(indices)}')

        plt.figure(figsize=(x, y))
        plt.pcolor(som.distance_map().T[:actual_grid_size[0], :actual_grid_size[1]], cmap='bone_r')  # Crop the displayed part
        plt.colorbar()
        plt.show()


gene_manager = GeneDataManager('localhost', 'root', '720027924', 'pineal_data')
selected_features = ['embryonic_day_21', 'postnatal_day_5', 'postnatal_day_20', 'postnatal_day_40',
                     'embryonic_night_21', 'postnatal_night_5', 'postnatal_night_20', 'postnatal_night_40']

SOMClusterer.perform_som_clustering(gene_manager, selected_features)



