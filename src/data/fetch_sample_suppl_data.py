""" Script to fetch supplemental informations about all the samples that are linked with the geo series of interest"""
import numpy as n
import pandas as pd
import pickle
import os
import sys
import re
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
import time
try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError       # for Python 2

print("\n## fetch detailed informations about all samples that are contained in the already fetched geo series ##\n")

## Set paths
# cdir = dir of this script
cdir = os.path.dirname(os.path.realpath(__file__))
# basedir = root dir of the repository
basedir = os.path.dirname(os.path.dirname(cdir))


#os.path.getmtime(path)


dir_data_in = basedir+"/data/interim/records_samples/"
dir_data_out = dir_data_in+"samples_suppl/"

# filenamebase for the large file containing all possibly relevant supplemental informations about the samples
fnameb_suppl = """samples_suppl"""
# filenamebase for the smaller file containing the informations about the sample title and the characteristics of the sample channels
fnameb_simple = 'samples_suppl_simple'

# if files were already writtn before, this variable is overwritten with the last geo series accession number from which samples have been read
last_gse_acc = ''
# number of already written files
fN = -1

if not os.path.exists(dir_data_out):
    os.makedirs(dir_data_out)
else:
    # check whether files have been produced before
    pattern = re.compile(r"""^"""+fnameb_suppl+"""_(?P<N>\d*).pkl$""")
    fN_list = []
    for f in os.listdir(dir_data_out):
        match = pattern.match(f)
        if(match):
            fN_list.append(int(match.group("N")))
    if len(fN_list) > 0:
        fN = max(fN_list)
        last_file = pickle.load( open(dir_data_out+fnameb_suppl+'_'+str(fN)+'.pkl','rb') )
        last_gse_acc = last_file['start_end_gse_acc'][-1]

def import_df(fname):
    print("start load samples database")
    df = pickle.load( open(dir_data_in+fname, 'rb') )
    print("done")
    return df

# import records dataframe
df_records = import_df('records.pkl')

if len(last_gse_acc) > 0:
    start_ndx = df_records.index[df_records['Accession'] == last_gse_acc].to_list()[0]
    if start_ndx != df_records.index[-1]:
        df_records = df_records.iloc[(start_ndx+1):,:]
    else:
        sys.exit("\nAll suppl. GEO samples data has been fetched and their relevant data been stored\n")


# url base from which we query
urlbase = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?'

# pattern to read the lines of the fetched data:
# tries to catch the key and the corresponding value of each line, both as strings
pattern = re.compile(r"""^[!\^]Sample[\_]*(?P<key>\D*?)(_ch)*(?P<ch>\d*)\s*=\s*(?P<value>.*)$""", re.IGNORECASE)

# the list of channel specific keys of samples for which we want to extract the values
ch_key_list = ['source_name', 'characteristics', 'treatment_protocol', 'growth_protocol', 'molecule', 'extract_protocol']

# the list of general keys of samples for which we want to extract the value
gn_key_list = ['title', 'status', 'submission_date', 'type', 'channel_count', 'description', 'data_processing']

# keys for which we try to flatten their list entries by joining them as strings
keys_flt_vals = ['treatment_protocol', 'extract_protocol', 'growth_protocol', 'data_processing']


# global list containing all sample accession numbers of all samples that will be fetched
glob_samples_acc_key_list = []

# global list containing all sample accession numbers of all samples that will be fetched. This one  will be reset after each batchsize number of samples have been fetched
glob_samples_acc_key_list_reset = []

# global list containing as elements all possibly relevant entries of the samples as dicts of all samples that will be fetched
glob_samples_dict_list = []

# global list containing as elements the name of the samples and the channel characteristics as dicts of all samples that will be fetched
glob_samples_dict_list_simple = []

# count the number of fetched samples, reset after each 10000 fetched samples
cnt_samples = 0

# global count the number of fetched samples, not reset
cnt_samples_glob = 0

# counter for files to be pickled
cnt_file = fN + 1


# the number of samples whose informations are contained in each part of the to be written samples_supp dict
batchsize = 10000

# max number of samples to be fetched, if <= 0, fetch all available samples
max_samples = -1

# flag to break the for loops when max_samples samples have been fetched
flag_break = False

# flag for determining if batchsize of samples have been reached, use -1 to start the loop, 0 if within the loop the batchsize has not been reached and 1 else
reached_batchsize = -1


class Sample:
    """ class to generate sample output from a list of strings that belong to the sample entry in the corresponding geo series and a regular expression pattern used for matching key, ch(annel number) and, (key) val(ue)
        - output is
            - self.out: dict of the gn_key_list keys and corresponding values together with the channel entries as further nested dicts
            - self.out_simple: dict with title as key and corresponding value together with numbers as additional keys identifying the channels whose values are the channel characteristics
    """
    def __init__(self, sample_lines, pattern):
        """ intit class """
        self.__lines = sample_lines[1:]
        self.__loc_sample_nest_key_list = gn_key_list[:]
        self.__loc_sample_nest_val_list = [[] for i in range(len(gn_key_list))]
        self.key = ''
        self.key_old = ''
        self.ch_old = -1
        self.ch = -1
        self.val = -1
        self.__parse_lines(pattern)
        self.__out()

    def __parse_lines(self, pattern):
        """ parse the lines in and match key, ch(annel number) and, (key) val(ue) with the help of the regular expression pattern """
        for line in self.__lines:
            match = pattern.match(line)
            if(bool(match)):
                self.key = match.group('key')
                self.ch  = match.group('ch')
                self.val = match.group('value')
                self.__clean()
                if len(self.ch) > 0:
                    self.ch = int(self.ch.strip())
                else:
                   self.ch = 0
                if self.key in gn_key_list:
                    if self.key != 'channel_count':
                        self.__loc_sample_nest_val_list[gn_key_list.index(self.key)].append(self.val)
                    else:
                        self.val = self.val.strip()
                        try:
                            self.val = int(self.val)
                        except:
                            pass
                        self.__loc_sample_nest_val_list[gn_key_list.index(self.key)] = self.val
                else:
                    if self.ch > 0:
                        if self.key in ch_key_list:
                            if self.ch not in self.__loc_sample_nest_key_list:
                                self.__loc_sample_nest_key_list.append(self.ch)
                                self.__loc_sample_nest_val_list.append([[] for i in range(len(ch_key_list))])
                                self.__loc_sample_nest_val_list[self.__loc_sample_nest_key_list.index(self.ch)][ch_key_list.index(self.key)].append(self.val)
                            else:
                                self.__loc_sample_nest_val_list[self.__loc_sample_nest_key_list.index(self.ch)][ch_key_list.index(self.key)].append(self.val)
                self.key_old = self.key
                self.ch_old = self.ch
            self.__clean()

    def __clean(self):
        """ clean the values list entries by joining string list entries to a single string list entry for certain keys """
        if( self.key_old != self.key ):
            if self.key_old in gn_key_list:
                if isinstance(self.__loc_sample_nest_val_list[gn_key_list.index(self.key_old)], list):
                    l_tmp_list = len(self.__loc_sample_nest_val_list[gn_key_list.index(self.key_old)])
                    if l_tmp_list > 1:
                        if self.key_old in keys_flt_vals:
                            self.__loc_sample_nest_val_list[gn_key_list.index(self.key_old)] = [" ".join(self.__loc_sample_nest_val_list[gn_key_list.index(self.key_old)])]
            else:
                if self.ch_old > 0:
                    if self.key_old in ch_key_list:
                        if self.ch_old in self.__loc_sample_nest_key_list:
                            if isinstance(self.__loc_sample_nest_val_list[self.__loc_sample_nest_key_list.index(self.ch_old)][ch_key_list.index(self.key_old)], list):
                                l_tmp_list = len(self.__loc_sample_nest_val_list[self.__loc_sample_nest_key_list.index(self.ch_old)][ch_key_list.index(self.key_old)])
                                if l_tmp_list > 1:
                                    if self.key_old in keys_flt_vals:
                                        self.__loc_sample_nest_val_list[self.__loc_sample_nest_key_list.index(self.ch_old)][ch_key_list.index(self.key_old)] = [" ".join(self.__loc_sample_nest_val_list[self.__loc_sample_nest_key_list.index(self.ch_old)][ch_key_list.index(self.key_old)])]

    def __out(self):
        """ method to store the relevant sample informations for writing """
        ch_count_list = list(set(self.__loc_sample_nest_key_list).difference(set(gn_key_list)))
        self.__loc_sample_nest_key_list_simple = ['title']
        self.__loc_sample_nest_val_list_simple = [self.__loc_sample_nest_val_list[self.__loc_sample_nest_key_list.index('title')]]
        if len(ch_count_list) > 0:
            for self.ch in ch_count_list:
                self.__loc_sample_nest_val_list[self.__loc_sample_nest_key_list.index(self.ch)] = dict(zip(ch_key_list, self.__loc_sample_nest_val_list[self.__loc_sample_nest_key_list.index(self.ch)]))
                self.__loc_sample_nest_key_list_simple.append(self.ch)
                self.__loc_sample_nest_val_list_simple.append(self.__loc_sample_nest_val_list[self.__loc_sample_nest_key_list.index(self.ch)]['characteristics'])
        self.out = dict(zip(self.__loc_sample_nest_key_list, self.__loc_sample_nest_val_list))
        self.out_simple = dict(zip(self.__loc_sample_nest_key_list_simple, self.__loc_sample_nest_val_list_simple))



for acc_series in df_records['Accession']:
    print(acc_series)
    # urllib request string: acc_series is the geo series acc number
    #   tries to fetch all informations in the SOFT(form=text) format
    request = 'acc='+acc_series+'&targ=all&form=text&view=brief'

    if(reached_batchsize in [-1,1]):
        start_acc_series = acc_series
        reached_batchsize = 0
    code=404
    while(code == 404 or code == 501):
        try:
            urllib_query = urlopen(urlbase, request.encode())
        except HTTPError as e:
            code = e.code
            print("HTTPError: "+str(code))
            time.sleep(1)
        except URLError as e:
            code = e.code
            print("URLError: "+str(code))
            time.sleep(1)
        else:
            code = urllib_query.getcode()
    gse_str_entries_raw = urllib_query.read()
    urllib_query.close()
    gse_str_entries_raw = gse_str_entries_raw.decode('utf-8', 'replace' )
    gse_str_entries = gse_str_entries_raw.splitlines()
    ndxl = [i for i, x in enumerate(gse_str_entries) if x.startswith('^SAMPLE')]
    ndxl.append(len(gse_str_entries))
    for i in range(len(ndxl[0:-1])):
        sample_lines = gse_str_entries[ndxl[i]:ndxl[i+1]]
        match = pattern.match(sample_lines[0])
        sample_acc_id = match.group('value')
        if sample_acc_id not in glob_samples_acc_key_list:
            glob_samples_acc_key_list.append(sample_acc_id)
            glob_samples_acc_key_list_reset.append(sample_acc_id)
            cnt_samples += 1
            cnt_samples_glob += 1
            print((cnt_samples_glob,sample_acc_id))

            sample = Sample(sample_lines, pattern)

            glob_samples_dict_list.append(sample.out)
            glob_samples_dict_list_simple.append(sample.out_simple)

            # break from inner loop when max_samples have been fetched
            if max_samples > 0:
                if cnt_samples_glob == max_samples:
                    flag_break = True
                    break

        # break from loop over sample text paragraphs, if flag_break is true
        if flag_break:
            break

    if cnt_samples > batchsize:
        # full filename for the part of the sample_supp dictionary to be stored to a file
        glob_samples_acc_key_list_reset.append('start_end_gse_acc')
        glob_samples_dict_list.append([start_acc_series, acc_series])
        glob_samples_dict_list_simple.append([start_acc_series, acc_series])
        fname_suppl = fnameb_suppl+'_'+str(cnt_file)+'.pkl'
        # write to pickle file part of the detailed info sample dict
        f = open(dir_data_out+fname_suppl,"wb")
        pickle.dump(dict(zip(glob_samples_acc_key_list_reset, glob_samples_dict_list)),f)
        f.close()
        # write to pickle file part of the simple info sample dict
        fname_simple = fnameb_simple+'_'+str(cnt_file)+'.pkl'
        # write to pickle file
        f = open(dir_data_out+fname_simple,"wb")
        pickle.dump(dict(zip(glob_samples_acc_key_list_reset, glob_samples_dict_list_simple)),f)
        f.close()
        # reset the global lists of acc keys and possibly relevant entries for all fetched samples so far
        glob_samples_acc_key_list_reset = []
        glob_samples_dict_list = []
        glob_samples_dict_list_simple = []
        # reset counter for samples
        cnt_samples = 0
        # increase counter that indicates the part of the samples_supp dict to be written to a file
        cnt_file += 1
        reached_batchsize = 1

    # break from loop over geo series, if flag_break is true
    if flag_break:
        break

if 'start_end_gse_acc' not in glob_samples_acc_key_list_reset:
    glob_samples_acc_key_list_reset.append('start_end_gse_acc')
    glob_samples_dict_list.append([start_acc_series, acc_series])
    glob_samples_dict_list_simple.append([start_acc_series, acc_series])

# write last part of the sample_supp dictionary to a file
fname_suppl = fnameb_suppl+'_'+str(cnt_file)+'.pkl'
# write to pickle file
f = open(dir_data_out+fname_suppl,"wb")
pickle.dump(dict(zip(glob_samples_acc_key_list_reset, glob_samples_dict_list)),f)
f.close()

# write the sample_supp dictionary to a file
fname_simple = fnameb_simple+'_'+str(cnt_file)+'.pkl'
# write to pickle file
f = open(dir_data_out+fname_simple,"wb")
pickle.dump(dict(zip(glob_samples_acc_key_list_reset, glob_samples_dict_list_simple)),f)
f.close()
