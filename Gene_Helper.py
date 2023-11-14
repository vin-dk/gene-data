import mysql.connector
import pandas as pd
import numpy as np 


# this is a helper class. It imports data and handles base handling operations. Provided the data follows a similar format, the variables can be changed to match up with your SQL database

class GeneData:
    def __init__(self, gene, gene_name, gene_biotype, chrom, start, stop,
                 embryonic_day_21, embryonic_night_21, postnatal_day_5, postnatal_night_5,
                 postnatal_day_20, postnatal_night_20, postnatal_day_40, postnatal_night_40,
                 mean, stddev, official_full_name, summary):
        self.gene = gene
        self.gene_name = gene_name
        self.gene_biotype = gene_biotype
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.embryonic_day_21 = embryonic_day_21
        self.embryonic_night_21 = embryonic_night_21
        self.postnatal_day_5 = postnatal_day_5
        self.postnatal_night_5 = postnatal_night_5
        self.postnatal_day_20 = postnatal_day_20
        self.postnatal_night_20 = postnatal_night_20
        self.postnatal_day_40 = postnatal_day_40
        self.postnatal_night_40 = postnatal_night_40
        self.mean = mean
        self.stddev = stddev     
        self.official_full_name = official_full_name
        self.summary = summary

class GeneDataManager:
    def __init__(self, host, user, password, database):
        self.host = host
        self.user = user
        self.password = password
        self.database = database

    def fetch_gene_data(self):
        conn = mysql.connector.connect(
            host=self.host,
            user=self.user,
            password=self.password,
            database=self.database
        )
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM pineal_data")
        gene_data_list = [GeneData(*row) for row in cursor]
        cursor.close()
        conn.close()
        return gene_data_list
    
    def gene_data_manager_merged(self):
        gene_data_list = self.fetch_gene_data()

        x_values = [getattr(gene_data, 'postnatal_day_40') for gene_data in gene_data_list]
        y_values = [getattr(gene_data, 'postnatal_day_20') for gene_data in gene_data_list]

        group1_avg = self.calculate_average(x_values)
        group2_avg = self.calculate_average(y_values)

        log_x_values = [self.log_transformation(value) for value in x_values]
        log_y_values = [self.log_transformation(value) for value in y_values]

        data = pd.DataFrame({'x_value': log_x_values, 'y_value': log_y_values})

        return self.significant_expression(data)    
    
    def calculate_average(self, values):
        return sum(values) / len(values)
    
    def log_transformation(self, value):
            return np.log10(value + 1)
    
    def significant_expression(self, data, threshold=1):
            significant_expression = (
                ((data['x_value'] == 0) & (data['y_value'] >= 0.5)) |
                ((data['y_value'] == 0) & (data['x_value'] >= 0.5)) |
                ((data['x_value'] > 0) & (data['y_value'] > 0) &
                 ((data['x_value'] >= (1 + threshold) * data['y_value']) | 
                  (data['y_value'] >= (1 + threshold) * data['x_value'])))
            )
            significant_data = data[significant_expression]
            return significant_data
    
    def analyze_significant_data(self, significant_data, gene_data_list):
            for index, row in significant_data.iterrows():
                x_value = row['x_value']
                y_value = row['y_value']
                gene_name = gene_data_list[index].gene_name
                multiplier = max(x_value, y_value) / min(x_value, y_value)
    
                print(f"Gene Name: {gene_name}")
                print(f"X Value: {x_value}")
                print(f"Y Value: {y_value}")
                print(f"X and Y are {multiplier:.2f} times different")    

    @staticmethod
    def get_group_values(gene_data_list, group_periods):
        group_values = {period: [getattr(gene_data, period) for gene_data in gene_data_list] for period in group_periods}
        return group_values

