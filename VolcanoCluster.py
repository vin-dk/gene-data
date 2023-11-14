import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Gene_Helper import GeneDataManager  
from scipy import stats

class VolcanoCluster:
    def __init__(self):
        pass

    def generate_volcano_plot(self):
        gene_manager = GeneDataManager('localhost', 'root', '720027924', 'pineal_data')
        gene_data_list = gene_manager.fetch_gene_data()

        group1_periods = ['embryonic_day_21', 'postnatal_day_5', 'postnatal_day_20', 'postnatal_day_40']
        group2_periods = ['embryonic_night_21', 'postnatal_night_5', 'postnatal_night_20', 'postnatal_night_40']

        p_val_ttest_dict = {}
        pseudocount = 0.01
        altered_max = 0
        fc_max = 0
        fc_min = 100000000

        for gene_data in gene_data_list:
            group1_values = [getattr(gene_data, period) for period in group1_periods]
            group2_values = [getattr(gene_data, period) for period in group2_periods]

            group1_avg = (sum(group1_values)) / len(group1_values) + pseudocount  # Average of group 1
            group2_avg = (sum(group2_values)) / len(group2_values) + pseudocount  # Average of group 2
            fold_change = (group2_avg / group1_avg)
            fold_change = np.log2(fold_change)

            if fold_change > fc_max:
                fc_max = fold_change
            if fold_change < fc_min:
                fc_min = fold_change
    
            if all(val == 0 for val in group1_values) and all(val == 0 for val in group2_values):
                p_val_ttest = 0  # Set p-value to 0 if all values are 0 (handle exception)
            else:
                p_val_ttest = stats.ttest_rel(group1_values, group2_values).pvalue
                p_val_ttest = -np.log10(p_val_ttest)
                if p_val_ttest > altered_max:
                    altered_max = p_val_ttest
    
            p_val_ttest_dict[gene_data.gene] = p_val_ttest
    
            gene_data.fold_change = fold_change
            gene_data.p_val_ttest = p_val_ttest
    
        significance_threshold = 1.3  # -log10(0.05) = 1.3
    
        # Prepare the data for plotting
        data_calculated = pd.DataFrame([vars(gene_data) for gene_data in gene_data_list])
    
        column_order = ['gene', 'gene_name', 'gene_biotype', 'chrom', 'start', 'stop',
                        'embryonic_day_21', 'embryonic_night_21', 'postnatal_day_5', 'postnatal_night_5',
                        'postnatal_day_20', 'postnatal_night_20', 'postnatal_day_40', 'postnatal_night_40',
                        'mean', 'stddev', 'fold_change', 'p_val_ttest',
                        'official_full_name', 'summary']

    
        data_reordered = data_calculated[column_order]
    
        plt.figure(figsize=(10, 6))
        scatter = plt.scatter(data_reordered['fold_change'], data_reordered['p_val_ttest'], c='grey', label='Other')
    
        plt.axvline(-2, color='k', linestyle='--', label='x = -2')
        plt.axvline(2, color='k', linestyle='--', label='x = 2')
        plt.axhline(significance_threshold, color='r', linestyle='--', label=f'y = {significance_threshold}')
    
        for i, row in data_reordered.iterrows():
            x = row['fold_change']
            y = row['p_val_ttest']
    
            if y >= significance_threshold:
                if x <= -2 or x >= 2:
                    plt.scatter(x, y, c='darkblue', marker='o')
                    gene_label = row['gene_name']  # display gene name
                    plt.annotate(gene_label, (x, y), textcoords="offset points", xytext=(0, 10), ha='center', rotation=45)
                else:
                    plt.scatter(x, y, c='lightblue', marker='o')
    
        plt.axhline(0, color='k', linestyle='--')
        plt.xlabel('Fold Change (log2)')
        plt.xlim(-15, 15)
        plt.ylabel('-log10(p-value)')
        plt.title('Day vs Night Gene Expression')
        plt.legend()
        plt.show()
        
volcano_cluster = VolcanoCluster()
volcano_cluster.generate_volcano_plot()