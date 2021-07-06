import pandas as pd

original = pd.read_csv('./x50_10Hz.txt',sep='\t')
original.drop_duplicates(keep=False,inplace=True)
original.to_csv('./x50_10Hz_filtered.txt',sep='\t')