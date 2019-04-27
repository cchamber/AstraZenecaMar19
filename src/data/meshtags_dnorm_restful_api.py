import time
import sys
import getopt
import urllib
from urllib.error import HTTPError, URLError
from urllib.request import urlopen
import numpy as np
import os
import pandas as pd
import glob
from shutil import copyfile as copy
import subprocess
import string


print("\n## obtain mesh ids for geo series by calling the DNORM program via the RESTFUL API ##\n")
print("\n## \t Note: This can take up to a couple of hours for many geo series ##\n")

## Set paths
# cdir = dir of this script
cdir = os.path.dirname(os.path.realpath(__file__))
# basedir = root dir of the repository
basedir = os.path.dirname(os.path.dirname(cdir))



dir_data_in = basedir+'/data/interim/records_samples/'
dir_data_out = basedir+"/data/interim/"
tagPath = basedir+'/data/interim/Tags'

tagInputFile = 'input_restful.txt'
tagInputFile = os.path.join(tagPath, tagInputFile)


tagOutputFile = '_output_restful.txt'
trigger = ['DNorm', 'tmChem']

if not os.path.exists(dir_data_out):
    os.makedirs(dir_data_out)

if not os.path.exists(tagPath):
    os.makedirs(tagPath)

# make input from records
# this cell makes input for taggers

# This cell creates input in the form of
# 10192393|t|A common human skin tumour is caused by activating mutations in beta-catenin.
# 10192393|a|WNT signalling orchestrates a number of developmental programs. 

df_record = pd.read_pickle(os.path.join(dir_data_in,'records.pkl'))

all_rows = 1
n_rows = 10

if all_rows==1:
    n_rows = len(df_record)

rec4input = df_record.loc[:n_rows,:]

outF = open(tagInputFile, "w")
# build list of strings
rec4input = rec4input.reset_index()
for i in range(len(rec4input)):
    Id = rec4input.loc[i,:].Id
    title = rec4input.loc[i,:].title.translate(str.maketrans('', '', string.punctuation))
    summary = rec4input.loc[i,:].summary.translate(str.maketrans('', '', string.punctuation))
    # remove | characters from title and summary
    outF.write(Id+'|t|'+title)
    outF.write("\n")
    outF.write(Id+'|a|'+summary)
    outF.write("\n")

# tag text summaries

print('This script uses NCBI restful API https://www.ncbi.nlm.nih.gov/research/bionlp/APIs/')
print('First we request tags from records in batches (e.g. n=1000), then we save the tags.\
      \n~45000 summaries should take ~30 minutes to tag with one tagger.')

def submit_request(url_Submit, InputSTR):
    urllib_submit = urllib.request.urlopen(url_Submit, InputSTR.encode())
    urllib_result = urllib.request.urlopen(url_Submit, InputSTR.encode())
    SessionNumber = urllib_submit.read()
    SessionNumber = SessionNumber.decode('utf-8')
    print("Thanks for your submission. The session number is : "+ str(SessionNumber))
    print("The request is received and processing....\n\n")
    code = urllib_result.getcode()
    return SessionNumber, code

def main_function(inputfile, outputfile, trigger):

    taxonomy = ''
    email = ''
    PubTator_username = ''
    url_Submit = ''

    if taxonomy != '':
        url_Submit = "https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/" + trigger + "/" + taxonomy + "/"
    elif email != '':
        url_Submit = "https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/" + trigger + "/Submit:" + email + "/"
    else:
        url_Submit = "https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/" + trigger + "/Submit/"

    fh = open(inputfile)
    InputLIST=[]

    for line in fh:    
        if line[9:12] == '|t|':
            zrec = line
        if line[9:12] == '|a|': 
            zrec = zrec+line
            InputLIST += [zrec]

    n_records = len(InputLIST)
    print('There are '+str(n_records)+' records.\n')
    batch_size = 2000
    n_groups = np.ceil(n_records/batch_size).astype(int)
    start = np.arange(0,n_records,batch_size).astype(int)

    tic = time.time() 
    left_over_list = []
    max_attempts = 2

    SessionList = []
    for i in np.arange(0,n_groups,1):#range(n_groups):
        print('Requesting tags for records from '+str(start[i])+' to '+str(np.min([n_records,start[i]+batch_size])))
        InputSTR = ''.join(InputLIST[start[i]:np.min([n_records,start[i]+batch_size])])
        code = 403
        n_attempt = 0
        while code == 403:
            try:
                SessionNumber,code = submit_request(url_Submit,InputSTR)
                SessionList +=[SessionNumber]
            except HTTPError as e:
                n_attempt += 1
                code = e.code
                print('http error')
                print(code)
                if n_attempt>max_attempts:
                    code = 1
                    left_over_list += [InputSTR]
            time.sleep(2)


    F = open(outputfile,'w')
    for iSessionNumber in SessionList:
        url_Receive = "https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/" + str(iSessionNumber) + "/Receive/"
        print('Attempt to save this URL to file (Do not click): '+url_Receive)
        code1=404
        while(code1 == 404 or code1 == 501):
            try:
                urllib_result = urllib.request.urlopen(url_Receive)
            except HTTPError as e:
                code1 = e.code
            except URLError as e:
                code1 = e.code
            else:
                code1 = urllib_result.getcode()
            toc = time.time()-tic
            print('Time elapsed since tagging began is ',str(np.round(toc*100)/100))
            print('code='+str(code1))
            time.sleep(5)  
        F.write(urllib_result.read().decode('utf-8')) 
    F.close()
    

for itrigger in trigger:
    print('\n\n')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('| Tagging with '+itrigger)
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
    main_function(tagInputFile, os.path.join(tagPath,itrigger+tagOutputFile), itrigger)

# save text files as pandas df

zdf = pd.DataFrame()

for itrigger in trigger:
    file_to_load = os.path.join(tagPath,itrigger+tagOutputFile)
    with open(file_to_load,'r') as f:
        data = f.read().splitlines()

    table = [d.split('\t') for d in data if len(d.split('\t'))>1]
    headers = ['Id', 'start', 'end', 'disease_tag', 'tag_type', 'ONTid']
    df = pd.DataFrame(table, columns=headers)

    s = df.ONTid.str.split(':')
    df_aux = pd.DataFrame.from_dict(dict(zip(s.index, s.values))).T
    df_aux.columns = ['ont', 'unique_id']
    df['ont'] = df_aux['ont']
    df['unique_id'] = df_aux['unique_id']
    df = df.drop('ONTid', axis = 1)
    zdf = pd.concat([zdf,df], axis=0)
zdf = zdf[zdf.ont!='None']
# delete mesh headings
zdf.to_pickle(os.path.join(dir_data_out,'geoid_meshui_dnorm_restful_api.pkl'))
