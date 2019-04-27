#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: nLp ATTACK

This code downloads tree IDs and supplemental concept records (SCRs) and puts them into one df

output df columns are: 
'category': record type (eg disease='C')
'mesh_id': unique id (UI) given to record (eg, D013614)
'mesh_heading': concept name (eg, 'Adams Nance syndrome')
'mesh_treenumbers': address in tree, with sub categories separated by . 
(eg C14.280.067.845.695), one concept can have many mesh_treenumbers
'scr': is this a supplemental concept record. 0 for no, 1 for yes. There are ~600K SCRs and ~50K MeSH records


SCRs Supplemental concept records:
more info: search MeSH here https://meshb.nlm.nih.gov/search and more info here https://www.nlm.nih.gov/mesh/intro_record_types.html
download here: ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/
Fields of SCR text file:
The names (NM) are new concepts
the NMs map to certain concepts that are already in the MeSH (HM)
The SCRs have their own unique identifiers UIs as well

How can we deal with SCR?
These new concepts might be tags of Karsten's PMID data for example
We can map the new record names NMs to existing MeSH ids (UIs and tree IDs)
-- this is the only way of linking this data with the tree. Then we search SCRs as if they were mesh headings
result: in our use case, we can search the tags for disease names, including SCR names
"""


print("\n## download mesh ids and its corresponding mesh headings and put it into the BASE/data/final/mesh.pkl dataframe ##\n")

import pandas as pd
import pickle
import numpy as np
import urllib.request
import os

# define 
include_SCR = 1

if include_SCR == 0:
    output_file = 'mesh.pkl'
elif include_SCR == 1:
    output_file = 'mesh_wSCR.pkl'

# cdir = dir of this script
cdir = os.path.dirname(os.path.realpath(__file__))
# basedir = root dir of the repository
basedir = os.path.dirname(os.path.dirname(cdir))

dir_data_in = basedir+"/data/external/"
dir_data_out = basedir+"/data/final/"

if not os.path.exists(dir_data_out):
    os.makedirs(dir_data_out)

if not os.path.exists(dir_data_in):
    os.makedirs(dir_data_in)

def fetch_mesh_remotely():
    files = ['d2019.bin', 'c2019.bin']
    for f in files:
        try:
            fh = open(dir_data_in+f,'r')
            fh.close
        except FileNotFoundError:
            url = 'ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/asciimesh/'+f
            urllib.request.urlretrieve(url, dir_data_in+f)


fetch_mesh_remotely()

# initialize
tree_value = []
name_list = []
tree_number_list = []
id_list = []

with open(dir_data_in+'d2019.bin') as f:
    for line in f: # cycle through each line
        
        if line.startswith('MH = '): # name
            name_list.append(line[5:-1])
        
        if line.startswith('MN = '): # tree numbers
            # collect tree number for each line            
            tree_value_temp = line[5:-1] 
            # include last char \n because it will help to search each level of the tree
            # collect all tree numbers
            tree_value.append(tree_value_temp)
               
        if line.startswith('UI = '): # unique id
            tree_number_list.append(tree_value)
            tree_value = [] # initialize since all tree numbers are obtained
            id_list.append(line[5:-1])
        
df = pd.DataFrame.from_dict({'mesh_id':pd.Series(id_list),'mesh_heading':pd.Series(name_list), 'mesh_treenumbers':pd.Series(tree_number_list)})     

no_tree_number_df = df[df.mesh_treenumbers.str.len()==0]

no_tree_number_df.loc[:, 'mesh_treenumbers'] = pd.Series(np.array(len(no_tree_number_df)*[np.nan]), index=no_tree_number_df.index)

no_tree_number_df.loc[:, 'category'] = pd.Series(len(no_tree_number_df)*['Sex'], index=no_tree_number_df.index)

# expand list into columns
tags = df.mesh_treenumbers.apply(pd.Series)
cols = ['tag'+str(icol) for icol in tags.columns]
tags.columns = cols
tags['mesh_id'] = df.mesh_id
df = pd.merge(df,tags, on='mesh_id', how='inner')
# melt
df = pd.melt(df, id_vars = ['mesh_id','mesh_heading'], value_vars=cols)
df = df.drop('variable',axis=1)
df.columns = ['mesh_id','mesh_heading', 'mesh_treenumbers']
df['category'] = df.mesh_treenumbers.str[:1]
df = df.dropna()

if include_SCR ==1: 
    new_name_list = []
    maps_to = []
    maps_to_list = []
    id_list = []

    with open(dir_data_in+'c2019.bin') as f:
        for line in f: # cycle through each line

            if line.startswith('NM = '): # new name
                new_name_list.append(line[5:-1])

            if line.startswith('HM = '): # maps to names in MeSH           
                maps_to_temp = line[5:-1].split('/')[0].replace('*','')
                maps_to.append(maps_to_temp)

            if line.startswith('UI = '): # unique id
                maps_to_list.append(maps_to)
                maps_to = [] # initialize since all tree numbers are obtained
                id_list.append(line[5:-1])
    scr_df = pd.DataFrame.from_dict({'SCR_id':pd.Series(id_list),'SCR_heading':pd.Series(new_name_list),'mesh_heading':pd.Series(maps_to_list)})     
    # expand list into columns, then melt
    tags = scr_df.mesh_heading.apply(pd.Series)
    cols = ['tag'+str(icol) for icol in tags.columns]
    tags.columns = cols
    tags['SCR_id'] = scr_df.SCR_id
    scr_df = pd.merge(scr_df,tags, on='SCR_id', how='inner')
    scr_df = pd.melt(scr_df, id_vars = ['SCR_id','SCR_heading'], value_vars=cols)
    scr_df = scr_df.drop('variable',axis=1)
    scr_df.columns = ['SCR_id','SCR_heading', 'mesh_heading']
    scr_df = scr_df.dropna()

    # Let's map SCR info to the main df info
    # SCR_heading becomes new mesh_heading
    # set variable map_scrUI_to_meshUI__OR__use_scrUI_as_UI: 
    # we can map the scr UIs to the existing mesh UIs (0) 
    # or use the scr UIs as UI (1)
    map_scrUI_to_meshUI__OR__use_scrUI_as_UI =1

    if map_scrUI_to_meshUI__OR__use_scrUI_as_UI ==0: 
        scr_df = pd.merge(scr_df[['SCR_heading', 'mesh_heading']], df, on='mesh_heading', how='inner')
        scr_df = scr_df.drop('mesh_heading', axis = 1)
        scr_df = scr_df.rename(index=str, columns = {'SCR_heading':'mesh_heading'})
    elif map_scrUI_to_meshUI__OR__use_scrUI_as_UI ==1:
        scr_df = pd.merge(scr_df[['SCR_heading', 'mesh_heading', 'SCR_id']], df, on='mesh_heading', how='inner')
        scr_df = scr_df.drop(['mesh_heading','mesh_id'], axis = 1)
        scr_df = scr_df.rename(index=str, columns = {'SCR_heading':'mesh_heading', 'SCR_id':'mesh_id'})
    
    scr_df['scr'] = 1

    df['scr'] = 0

    df = pd.concat([df, scr_df], ignore_index=True, sort = False)

df = pd.concat([df, no_tree_number_df], ignore_index=True, sort = False)

df.to_pickle(os.path.join(dir_data_out,output_file))
