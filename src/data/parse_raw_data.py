""" XML parsing
In this file we will parse the xml data files

1. First we loop through xml files and get a very long list of dictionaries.
2. We turn the list of dictionaries into a python data frame.
3. Then we get the data for the sample, stored in subdictionaries in the 'Samples' field of the main dataframe.
4. Number of samples to the records df
5. Output is saved two dataframes, converted to pkl files, in the same folder as the raw data files:
    - df_records --> records.pkl
    - df_samples --> samples.pkl
"""
from Bio import Entrez
import pandas as pd
import numpy as np
import glob
import os
"""Parsing raw data
- First we loop through xml files and get a very long list of dictionaries.
- We turn the list of dictionaries into a pandas data frame.
- Then we go through each row of the larger dataframe and get the sample data from each row (takes a while)
"""
Entrez.email = "A.N.Other@example.com" # Always tell NCBI who you are

print("\n## parse data into dataframes ##\n")

## Set paths
# cdir = dir of this script
cdir = os.path.dirname(os.path.realpath(__file__))
# basedir = root dir of the repository
basedir = os.path.dirname(os.path.dirname(cdir))


dir_data_in = basedir+'/data/raw/geo_data'              # dir for reading
dir_data_out = basedir+'/data/interim/records_samples'  # dir for storing

if not os.path.exists(dir_data_out):
    os.makedirs(dir_data_out)

output_file_rec = os.path.join(dir_data_out, 'records.pkl')
output_file_sam = os.path.join(dir_data_out, 'samples.pkl')
raw_files = sorted(glob.glob(os.path.join(dir_data_in, '*.xml')))

if os.path.isdir(dir_data_out)!=1:
    os.mkdir(dir_data_out)

record_list = []
for ifile in raw_files:
    print('Parsing ',ifile)
    handle = open(ifile)
    records = Entrez.parse(handle)
    for record in records:
        record_list.append(record)

df_records = pd.DataFrame(record_list)

# generating the samples df
list_geoid_global = []
list_nsamples_global = []
list_date_global = []
list_accession_global = []
list_title_global = []

for i in df_records.index:
    list_accession_local = []
    list_title_local = []
    
    geoid = str(df_records.iloc[i].Id)
    n_samples = df_records.iloc[i].n_samples
    date = str(df_records.iloc[i].PDAT)

    for sample in df_records.loc[i].Samples:
        list_accession_local.extend([sample['Accession']])
        list_title_local.extend([sample['Title']])
    
    list_geoid_global.extend(n_samples*[geoid])
    list_nsamples_global.extend(n_samples*[n_samples])
    list_date_global.extend(n_samples*[date])

    list_accession_global.extend(list_accession_local)
    list_title_global.extend(list_title_local)

    if i%5000==0:
        print('Sample iteration:')
        print(i)

df_samples = pd.DataFrame({'geo_id': np.array(list_geoid_global), 'nsamples': np.array(list_nsamples_global), 'date': np.array(list_date_global), 'accession': np.array(list_accession_global), 'title': np.array(list_title_global)})

df_samples = df_samples[['geo_id', 'nsamples', 'date', 'accession', 'title']]

df_samples = df_samples.sort_values(by=['geo_id']).reset_index(drop = True)

print('Saving samples to ', output_file_sam)
df_samples.to_pickle(output_file_sam)

# delete the samples columns since this is now included in the samples df
df_records = df_records.drop(columns=['Samples'])
# process the records df for easier accessibility
df_records['Id'] = df_records['Id'].apply(lambda x: str(x))
df_records['PubMedIds'] = df_records['PubMedIds'].apply(lambda x: np.array(list(x)).astype(int))
df_records['PDAT'] = df_records['PDAT'].apply(lambda x: str(x))

print('')
print('Saving records to ', output_file_rec)
print('')
df_records.to_pickle(output_file_rec)

print('Done.')
