#### script with common functions that I keep using (hopefully I'll keep this updated :)


import numpy as np

def make_equal_hist_bins(indata, bins=10):
    # for different sized datasets, make the binning = in a histogram
    bin_edges = np.linspace(indata.min(), indata.max(), bins + 1)
    return bin_edges



#### for importing from text files (also pandas can do this automatically I think)
def get_ints(data,sepa=' ',reshape=(-1,2)):
    ''' Function to get integer or float values from excel data (usually tab or ' ' separated).'''
    return np.fromstring(data,sep=sepa).reshape(reshape)

def get_strings(data, sepa='\t'):
    ''' Function to get string values from excel data (usually tab or ' ' separated).'''
    return [i.split(sepa,1)[0] for i in data.strip().split('\n')]
    