
import numpy as np
import pandas as pd

import os
import time


working_directory = "/home/robbler/research/Analysis_graphs/"
details_filename = "table_4_GNGS_edited - tsv_export.tsv"
detailed_file = working_directory+details_filename


# get the last time the file was updated
ti_m = os.path.getmtime(detailed_file)
m_ti = time.ctime(ti_m)

def disclaimer():
    print(f"[!!!] Reference file ({details_filename}) was last updated at: {m_ti}")


table = pd.read_csv(detailed_file,sep='\t')
table.sort_values(by='ID',ascending=True,inplace=True)
pd.set_option('display.precision', 10)
pd.set_option('display.max_colwidth', None)

######## Script to import data because I keep copy pasting these lines
######## Meant to import from TSV, CSV, or ' 'SV data like what comes from auto copy-pasting columns from spreadsheets...
def get_ints(data,sepa=' ',reshape=(-1,2)):
    ''' Function to get integer or float values from excel data (usually tab or ' ' separated).'''
    return np.fromstring(data,sep=sepa).reshape(reshape)

def get_strings(data, sepa='\t'):
    ''' Function to get string values from excel data (usually tab or ' ' separated).'''
    return [i.split(sepa,1)[0] for i in data.strip().split('\n')]

def add_obs_filters(tab,search_z='JADES z'):
    
    ### now add the restframe wavelength for each given the proper filter for them ###
    filters_dict = {
    "f444w": (3.0, 4.0),  # (min, max)
    "f356w": (2.4, 3.0),
    "f277w": (1.5, 2.4),
    "f200w": (0.75, 1.5),
    "f150w": (0.0, 0.75)
    } 
    wv_dict = { # in micrometers from https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters#gsc.tab=0
    "f444w": 4.401,
    "f356w": 3.566,
    "f277w": 2.776,
    "f200w": 1.988,
    "f150w": 1.501
    }  
    tab['Obs Filter'] = "None"
    tab['Obs WV (um)'] = "None"
    for f in filters_dict:
        # np.where((tab['Mostly Disk (threshold)'] == "no") & (tab['Mostly Spheroid']== "no") & (tab['Mostly Irregular']== "no"),labels[7],tab['Classification'])
        tab['Obs Filter'] = np.where((tab[search_z]>filters_dict[f][0]) & (tab[search_z]<=filters_dict[f][1]),f,tab['Obs Filter'])
        tab['Obs WV (um)'] = np.where((tab[search_z]>filters_dict[f][0]) & (tab[search_z]<=filters_dict[f][1]),wv_dict[f],tab['Obs WV (um)'])

    #### now find the rest frame wavelength for each source
    tab['Rest WV (um)'] = tab['Obs WV (um)']/(1+tab[search_z])

    return tab

def get_visual_morphology():
    '''
    function to get results of visual morphology survey
    '''
    ### fill out this function later if needed...
    return

def get_sources(positions=True,agn=False,full=False):
    '''
    function to get all GN & GS sources (names & locations)
    returns (kirk name, jades id, position)
    '''
    disclaimer()

    add_obs_filters(table) # assuming the data has 'JADES z' for the z column searching
    if(full):
        table_reduced=table
    elif(agn):
        table_reduced = table[['ID','JADES ID','RA Adjusted (Kirk)','Decl. Adjusted (Kirk)','AGN^a','z^b']]
    elif(positions):
        table_reduced = table[['ID','JADES ID','RA Adjusted (Kirk)','Decl. Adjusted (Kirk)']]
    else:
        table_reduced = table[['ID','JADES ID']]
    return table_reduced

def get_GN():
    ''' function to get GN names and locations
        returns (kirk name, jades id, location)
    '''

    return get_sources(full=True)

def get_GS():
    ''' function to get GS names and locations
        returns (kirk name, jades id, location)
    '''
    
    return get_sources(full=True)