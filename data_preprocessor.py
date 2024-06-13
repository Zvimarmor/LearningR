import pandas as pd
import numpy as np


# # data = pd.read_csv('GSE208783_series_matrix.txt', sep='\t', skiprows=40)
# data = data.dropna()
# data = data.transpose()

# data = data.replace('!', '', regex=True)
# data = data.dropna(axis=1, how='all')

# data.columns = data.iloc[0]
# data = data[1:]
# data = data.reset_index(drop=True)

# print(data.columns)

# data.drop('Sample_data_processing', axis=1, inplace=True)
# data.drop('Sample_relation', axis=1, inplace=True)


# #write the data to a new txt file
# data.to_csv('processed_data1.txt', sep='\t', index=False)

# # Load the data from both files into pandas DataFrames
# gene_data1_df = pd.read_csv('gene_data1.txt', sep='\t')
# gene_data2_df = pd.read_csv('gene_data2.txt', sep='\t')

# # Rename 'Gene symbol' column to 'Sample_ID' in gene_data1_df
# gene_data1_df.rename(columns={'Gene symbol': 'Sample_ID'}, inplace=True)

# # Merge the two DataFrames on the 'Sample_ID' column
# combined_df = pd.merge(gene_data1_df, gene_data2_df, on='Sample_ID', how='outer')

# # Save the combined DataFrame to a new file
# combined_df.to_csv('combined_gene_data_matrix.txt', sep='\t', index=False)


data = pd.read_csv('gene_data_combined.txt', sep='\t')
data = data.dropna()
data.to_csv('gene_data_combined.txt', sep='\t', index=False)




