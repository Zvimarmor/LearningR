import pandas as pd
import numpy as np


data = pd.read_csv('GSE208783_series_matrix.txt', sep='\t', skiprows=40)
data = data.dropna()
data = data.transpose()

data = data.replace('!', '', regex=True)
data = data.dropna(axis=1, how='all')

data.columns = data.iloc[0]
data = data[1:]
data = data.reset_index(drop=True)

print(data.columns)

data.drop('Sample_data_processing', axis=1, inplace=True)
data.drop('Sample_relation', axis=1, inplace=True)


#write the data to a new txt file
data.to_csv('processed_data1.txt', sep='\t', index=False)






