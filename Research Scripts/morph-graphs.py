import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
import import_data


import os
import time


### we could automatically download the analysis, but it involves google api & cloud project blah blah too much time
### see disclaimer() fucntion
working_directory = "/home/robbler/research/Analysis_graphs/"
filename = "Galaxy Classifications (Responses) - CSV-Analysis.tsv"
details_filename = "table_4_GNGS_edited - tsv_export.tsv"
file = working_directory+filename
detailed_file = working_directory+details_filename



# get the last time the file was updated
# ti_m = os.path.getmtime(file)
# m_ti = time.ctime(ti_m)

def disclaimer(input_file):
    ti_m = os.path.getmtime(input_file)
    m_ti = time.ctime(ti_m)
    in_filename = input_file.split('/')[-1]
    print(f"[!!!] Reference file ({in_filename}) was last updated at: {m_ti}")

disclaimer(file)


##############################################
####### taken from Conselice 2014   ##########
####### this is for the optical R-band   #####
####### we're in the ~1 mircon NIR/ I-band? ##
conselice_values = {
# Galaxy type: [Concentration (R), <-- (+/-), Asymmetry (R), <-- (+/-), Clumpiness (R), <-- (+/-)]
"Ellipticals":[4.4,0.3,0.02,0.02,0.00,0.04],
"Early-type disks (Sa-Sb)":[3.9,0.5,0.07,0.04,0.08,0.08],
"Late-type disks (Sc-Sd)":[3.1,0.4,0.15,0.06,0.29,0.13],
"Irregulars":[2.9,0.3,0.17,0.10,0.40,0.20],
"Edge-on disks":[3.7,0.6,0.17,0.11,0.45,0.20],
"ULIRGs":[3.5,0.7,0.32,0.19,0.50,0.40],
"Starbursts":[2.7,0.2,0.53,0.22,0.74,0.25],
"Dwarf ellipticals":[2.5,0.3,0.02,0.03,0.00,0.06]
}


### set plot font sizes
small_size = 16
large_size = 18
matplotlib.rc('font', size=large_size)
matplotlib.rc('axes', titlesize=large_size)
matplotlib.rc('legend',fontsize=small_size)
# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def plot_CEERS_testing():
    ## plot the CEERS subsample from earlier (6-12) to see if the asymmetry vs clumpiness makes sense compare to our data
    ### scatter the asymmetry vs smoothness (clumpiness) measurements
    ceers_data = pd.read_csv('/home/robbler/research/CEERS_statmorph_testing/testing_new_CEERS_statmorph.csv')
    fig = plt.figure(figsize=(8,6))

    ## overplot the Conselice values in gray to show the different ranges
    # conselice_measurements.plot(x='Clumpiness (R)',y='Asymmetry (R)',kind='scatter',xerr='Clumpiness (+/-)',yerr='Asymmetry (+/-)',c='gray')
    for g in conselice_values:
        clumpi = conselice_values[g][4]
        clumpi_err = conselice_values[g][5]
        asymm = conselice_values[g][2]
        asymm_err = conselice_values[g][3]
        plt.errorbar(clumpi,asymm,xerr=clumpi_err,yerr=asymm_err,fmt='.',alpha=0.5,color='darkgray')
        plt.annotate(g,(clumpi,asymm))

    plt.scatter(ceers_data['smoothness (S)'],ceers_data['asymmetry (A)'])

    plt.xlabel('Clumpiness (S)')
    plt.ylabel('Asymmetry (A)')
    # plt.xlim([-.01,.25])
    # plt.legend(loc='best',prop={"size":8})
    plt.tight_layout(pad=2)

    plt.title(f'f200w CEERS A vs S')
    plt.gca().invert_yaxis()
    # plt.savefig(f'{working_directory}/Asymmetry_Clumpiness/A_S_scatter_with_conselice_full_{nir_filter}.png',dpi=300)
    plt.show()
    plt.close()

table = pd.read_csv(file,sep='\t')
# table.sort_values(by='ID',ascending=True,inplace=True)
table.sort_values(by='ID',ascending=True,inplace=True)



def remove_excluded(tab,name='ID'):
    #### Remove excluded sources from table ####
    ######## ideally they're already removed, but we might have some mismatched data, so removed just to be safe
    excluded_sources = np.array(["GS_IRS35","GS_IRS74","GS_IRS80^e","GS_IRS81","GN_IRS26","GN_IRS31^e"]) # & GS_IRS1 was removed before passing it to this
    for i in excluded_sources:
        ind = tab[(tab[name].str.startswith(i))].index
        tab.drop(ind,axis=0,inplace=True)
    print(f"[---] Removed sources from table: {name=}")

#remove_excluded(table)



### apparently pandas has a query method??
##### get the different major classifications (only one galaxy per classification type)
def add_classifications(tab):
    ### function to add a classifications column to the inputted dataframe based on having 'yes' or 'no' in each respective column
    disk_only = tab.query('`Mostly Disk (threshold)` == "yes" & `Mostly Spheroid` == "no" & `Mostly Irregular` == "no"')
    disk_irr = tab.query('`Mostly Disk (threshold)` == "yes" & `Mostly Spheroid` == "no" & `Mostly Irregular` == "yes"')
    disk_sph = tab.query('`Mostly Disk (threshold)` == "yes" & `Mostly Spheroid` == "yes" & `Mostly Irregular` == "no"')
    spheroid_only = tab.query('`Mostly Disk (threshold)` == "no" & `Mostly Spheroid` == "yes" & `Mostly Irregular` == "no"')
    irr_sph = tab.query('`Mostly Disk (threshold)` == "no" & `Mostly Spheroid` == "yes" & `Mostly Irregular` == "yes"')
    irregular_only = tab.query('`Mostly Disk (threshold)` == "no" & `Mostly Spheroid` == "no" & `Mostly Irregular` == "yes"')
    disk_irr_sph = tab.query('`Mostly Disk (threshold)` == "yes" & `Mostly Spheroid` == "yes" & `Mostly Irregular` == "yes"')
    orphans = tab.query('`Mostly Disk (threshold)` == "no" & `Mostly Spheroid` == "no" & `Mostly Irregular` == "no"')

    tab['Main Class'] = np.where((tab['Disk? (%)'] > tab['Spheroid? (%)']) & (tab['Disk? (%)'] > tab['Irregular? (%)']),'Disk','Even Agreement')
    tab['Main Class'] = np.where((tab['Spheroid? (%)'] > tab['Disk? (%)']) & (tab['Spheroid? (%)'] >= tab['Irregular? (%)']),'Spheroid',tab['Main Class']) 
    # >= for spheroid because GN_IRS7^d is rest fram ~1 micron in f277w and has even classification between spheroid and irregular, but appears spheroidal in f277w (irregular in f150w, but different light probing there ~0.5 micron)
    tab['Main Class'] = np.where((tab['Irregular? (%)'] > tab['Disk? (%)']) & (tab['Irregular? (%)'] > tab['Spheroid? (%)']),'Irregular',tab['Main Class'])

    labels = ['Mostly Disk', 'Mostly Spheroid', 'Mostly Irregular','Disk+Spheroid','Disk+Irregular','Irregular+Spheroid','Disk+Irregular+Spheroid','Unclassifiable']

    tab['Classification'] = np.where((tab['Mostly Disk (threshold)'] == "yes") & (tab['Mostly Spheroid']== "no") & (tab['Mostly Irregular']== "no"),labels[0],'Orphan')
    tab['Classification'] = np.where((tab['Mostly Disk (threshold)'] == "no") & (tab['Mostly Spheroid']== "yes") & (tab['Mostly Irregular']== "no"),labels[1],tab['Classification'])
    tab['Classification'] = np.where((tab['Mostly Disk (threshold)'] == "no") & (tab['Mostly Spheroid']== "no") & (tab['Mostly Irregular']== "yes"),labels[2],tab['Classification'])
    tab['Classification'] = np.where((tab['Mostly Disk (threshold)'] == "yes") & (tab['Mostly Spheroid']== "yes") & (tab['Mostly Irregular']== "no"),labels[3],tab['Classification'])
    tab['Classification'] = np.where((tab['Mostly Disk (threshold)'] == "yes") & (tab['Mostly Spheroid']== "no") & (tab['Mostly Irregular']== "yes"),labels[4],tab['Classification'])
    tab['Classification'] = np.where((tab['Mostly Disk (threshold)'] == "no") & (tab['Mostly Spheroid']== "yes") & (tab['Mostly Irregular']== "yes"),labels[5],tab['Classification'])
    tab['Classification'] = np.where((tab['Mostly Disk (threshold)'] == "yes") & (tab['Mostly Spheroid']== "yes") & (tab['Mostly Irregular']== "yes"),labels[6],tab['Classification'])
    tab['Classification'] = np.where((tab['Mostly Disk (threshold)'] == "no") & (tab['Mostly Spheroid']== "no") & (tab['Mostly Irregular']== "no"),labels[7],tab['Classification'])

    # ### add them to general array
    morph_types = np.array([disk_only,spheroid_only,irregular_only,disk_sph,disk_irr,irr_sph,disk_irr_sph,orphans],dtype=object)

    # classes = np.full(table['ID'].shape[0],fill_value='None',dtype=object)
    # z = 0
    # for i in morph_types:
    #     ind = i.index.tolist()
    #     classes[ind] = labels[z]
    #     z+=1

    # tab.sort_index(ascending=True,inplace=True)
    # tab['Classification'] = classes
    # tab.sort_values(by='ID',ascending=True,inplace=True)
    
    from import_data import add_obs_filters
    ### but this might be the wrong z values (from Kirkpatrick not JADES)
    ### i mean either way it shouldn't affect the distribution of the variance from of rest frame wavelength too much
    add_obs_filters(tab,search_z='z')
    # print(tab['z'])

    return morph_types,labels
morph_types, labels = add_classifications(table)


#### importing statmorph values here::
### import the statmorph measurements tsv files *********************************************************************************************************
stat_measures = f"research/statmorph_output/grizli/12-23-statmorph-output.tsv"
stats = pd.read_csv(stat_measures,sep='\t')
disclaimer(stat_measures)
stats.sort_values(by='ID',ascending=True,inplace=True)
table.sort_values(by='ID',ascending=True,inplace=True)

## combine visual with this nir_filter statmorph measurement data
# combined = pd.concat([table,stats],axis=1,join='outer')
combined = pd.merge(table, stats, on="ID")
combined.sort_values(by='ID',ascending=True,inplace=True)

### Remove all Flag>1 # some filters will have different results this way though, so maybe we need to manually clean them before input
keep = combined.query('`Flag`<=1')

# keep['Merger (flag threshold)'].map(dict(yes=1, no=0))
# keep['Merger (flag threshold)'].replace(('yes', 'no'), (1, 0), inplace=True)
combined = keep
# print(combined['Merger (flag threshold)'])

### discard the bimodal measurements
def remove_bad_measurements():
    ## put a nan on measurements that should be thrown out
    bad_meas = f"research/statmorph_output/grizli/errorbar-trouble-sources.tsv"
    bad_ones = pd.read_csv(bad_meas,sep='\t')
    disclaimer(bad_meas)
    print('Removed some measurements for the following sources (bc bimodal error distribution)...')
    for id in bad_ones['id']:
        # go to the measurement value column & set it =np.nan
        # measurement = bad_ones[bad_ones['id']==id]['measurement']
        # print(combined['Concentration (C)'])
        combined.loc[combined['ID']==id,bad_ones[bad_ones['id']==id]['measurement']]=np.nan
        # print('would have removed measurement, but not removing for AAS since no error bars (1/10/25)')
        # print(combined[combined['ID']==id][bad_ones[bad_ones['id']==id]['measurement']])
    print(f'{bad_ones[["id","measurement"]]}')

remove_bad_measurements()

##### get the merger variants for each classification
##### for now we're including the 'maybe' flag, but in the future we should exclude it (or just exclude if from the sheet)
# v = 0
mergers = np.empty(len(morph_types),dtype=object)
mergers = combined.query('`Merger (flag threshold)` == "yes"')
# mergers = table.query('`Merger (flag threshold)` == "yes"')
# for i in morph_types:
    # mergers[v] = i.query('`Merger (flag threshold)` == "yes" | `Merger (flag threshold)` == "maybe"')
    # v+=1


#### plotting arrays ##### unused as of 7/24
disk_agn = morph_types[0]['AGN(%)']
sph_agn = morph_types[1]['AGN(%)']
irr_agn = morph_types[2]['AGN(%)']

disk_sph_agn = morph_types[3]['AGN(%)']
disk_irr_agn = morph_types[4]['AGN(%)']
irr_sph_agn = morph_types[5]['AGN(%)']
disk_irr_sph_agn = morph_types[6]['AGN(%)']
orphans = morph_types[7]

print(f'Total # of sources that match criteria: {sum(i.shape[0] for i in morph_types)}')
# print(f'Total # of sources in original table: {table.shape[0]}')
print(f'Total # of sources in original table: {combined.shape[0]}')
print(f'Total orphan sources: {orphans.shape[0]} ({orphans["ID"]})')

#### show # in each category
c=0
for i in morph_types:
    print(f'{labels[c]}: {i.shape[0]}')
    c+=1

morph_agn = [i[['AGN(%)']] for i in morph_types]
morph_z = [i[['z']] for i in morph_types]
morph_clumps = [i.query('`Has Clumps (flag)`=="yes"') for i in morph_types]
morph_spiral = [i.query('`Spiral Arms (flag)`=="yes"') for i in morph_types]


### categories for the visual classification graphing
# symbols = ['*','o','^','D','x','h','P','s']
bin_colors = {
    'Mostly Disk': ['#1f77b4','*'],
    'Mostly Spheroid': ['#ff7f0e','o'],
    'Mostly Irregular': ['#2ca02c','^'],
    'Disk+Spheroid': ['#d62728','D'],
    'Disk+Irregular': ['#9467bd','x'],
    'Irregular+Spheroid': ['#8c564b','h'],
    'Disk+Irregular+Spheroid': ['#e377c2','P'],
    'Unclassifiable': ['#7f7f7f','s']
}

main_colors = {
    'Disk': ['#1f77b4','*'],
    'Spheroid': ['#ff7f0e','o'],
    'Irregular': ['#2ca02c','^'],
    'Even Agreement': ['#7f7f7f','s']
}

base_wv_rest = np.median(table['Rest WV (um)'])
wv_rest_16 = np.abs(base_wv_rest-np.percentile(table['Rest WV (um)'],16))
wv_rest_84 = np.abs(base_wv_rest-np.percentile(table['Rest WV (um)'],84))

# base_wv_rest = table['Rest WV (um)'].mean()
# err_wv_rest = table['Rest WV (um)'].std()

def rest_frame_plot():
    ##### plot each filter rest wavelength as a function of z
    ##### then plot each source using several rest filters??
    filter_options = [1.15,1.501,1.988,2.776,3.566,4.401]
    # redshift = np.arange(0,4.0,.2,dtype=float)
    redshift = np.linspace(0,4.0)

    labels = ['f115w','f150w','f200w','f277w','f356w','f444w']
    c=0

    plt.scatter(table['z'],filter_options[3]/(1+table['z']),marker='o')
    for i in filter_options:
    #    filter_rest[c] = i/(1+redshift)
       plt.plot(redshift,i/(1+redshift),label=labels[c])
    #    plt.scatter(table['z'],i/(1+table['z']),marker='o')
       c+=1
    plt.xlabel('Redshift')
    plt.ylabel('Rest $\lambda$ ($\mu$m)')
    # plt.title(f'Rest-Frame Matching ($\lambda_{{rest}}={base_wv_rest:.2f}\pm{err_wv_rest:.3f} \mu m)$')
    plt.title(f'Rest-Frame NIRCam Filter Matching')
    ## filter alignments
    # plt.axhline(.75,c='gray')
    # plt.axhline(base_wv_rest-err_wv_rest,c='gray',alpha=0.3)
    # plt.axhline(base_wv_rest,c='gray',alpha=0.4)
    # plt.axhline(base_wv_rest+err_wv_rest,c='gray',alpha=0.3)

    ## redshift bins (matching the 1.0 Rest Wv (um) line)
    # plt.axvspan(3.0,4.0,facecolor='brown',alpha=0.3) # f444w
    # plt.axvspan(2.4,3.0,facecolor='purple',alpha=0.3) # f356w
    # plt.axvspan(0,4.0,facecolor='red',alpha=0.3) # f277w (this ones kind of large and off)
    # plt.axvspan(0.75,1.5,facecolor='green',alpha=0.3) # f200w
    # plt.axvspan(0,0.75,facecolor='orange',alpha=0.3) # f150w
    
    # plt.annotate('f444w',(3.5,4.4),horizontalalignment='center',verticalalignment='top')
    # plt.annotate('f356w',(2.65,4.4),horizontalalignment='center',verticalalignment='top')
    # plt.annotate('f277w',(2.0,4.4),horizontalalignment='center',verticalalignment='top')
    # plt.annotate('f200w',(1.0,4.4),horizontalalignment='center',verticalalignment='top')
    # plt.annotate('f150w',(.4,4.4),horizontalalignment='center',verticalalignment='top')

    plt.legend(loc='center right')
    plt.grid(True)
    plt.xticks()
    plt.yticks()
    plt.savefig(f'{working_directory}/Rest_frame_testing_plot.png',dpi=200)
    plt.show()
def rest_wavelength_hist():
    ### plot histogram of the distribution of rest wavelengths (micrometers)
    plt.hist(table['Rest WV (um)'],bins=5)
    plt.vlines(base_wv_rest-wv_rest_16,0,15,color='gray',linestyles='dashed')
    plt.vlines(base_wv_rest,0,15,color='red',linestyles='dashed')
    plt.vlines(base_wv_rest+wv_rest_84,0,15,color='gray',linestyles='dashed')
    plt.title('$\lambda_{rest}$ Distribution ($\mu m$)')
    plt.xlabel('$\lambda_{rest} (\mu m)$')
    plt.annotate(f'$\lambda_{{rest}}={base_wv_rest:.2}^{{+{wv_rest_84:.2}}}_{{-{wv_rest_16:.2}}}$',(1.06,14))
    plt.savefig(f'{working_directory}/Rest_frame_distribution_histogram.png',dpi=200)
    plt.show()

# rest_frame_plot()
# rest_wavelength_hist()
# exit()

######## graphs ###########

def agn_frac_hist(): #### need to update this one if we're gonna keep using it ** (7/15)
    # maybe do a stacked hist with spirals and clumps for each classification?
    ##### Make a plot of galaxy type by agn fraction
    histbins = 5
    bins=np.histogram(np.hstack((disk_agn,sph_agn,irr_agn,disk_sph_agn,disk_irr_agn,disk_irr_sph_agn,orphans['AGN(%)'])), bins=histbins)[1] #get the bin edges

    bins = np.arange(0,110,10)
    labels = np.arange(0,100,10)  # Labels for the bins
    table['AGN (bin)'] = pd.cut(table['AGN(%)'], bins=bins, labels=labels, right=False)
# table.grouby
    agns = combined.groupby(['AGN (bin)','Classification']).size().unstack()
    agns.plot(kind='bar',stacked=True, figsize=(10,6))


    # plt.hist(disk_agn, label='Mostly Disk',bins=bins,ec='blue',lw=1,fill=False,histtype='step')
    # plt.hist(sph_agn, label='Mostly Spheroid',bins=bins,ec='red',lw=1,fill=False)
    # plt.hist(irr_agn, label='Mostly Irregular',bins=bins,ec='green',lw=1,fill=False)
    # plt.hist(disk_sph_agn, label='Disk+Spheroid',bins=bins,ec='magenta',lw=1,fill=False)
    # plt.hist(disk_irr_agn, label='Disk+Irregular',bins=bins,ec='orange',lw=1,fill=False)
    # plt.hist(irr_sph_agn, label='Irregular+Spheroid',bins=bins,ec='yellow',lw=1,fill=False)
    # plt.hist(disk_irr_sph_agn, label='Disk+Spheroid+Irregular',bins=bins,ec='cyan',lw=1,fill=False)
    # plt.hist(orphans['AGN(%)'],label='Unclassifiable',bins=bins,ec='gray',lw=1,fill=False)

    plt.title("AGN Fraction by Galaxy Type")
    plt.ylabel('Count')
    plt.xlabel('AGN Fraction (%)')
    plt.legend()
    # plt.vlines([20,80],0,8)
    agn_filename = 'agn_frac_hist_all'
    plt.savefig(f'{working_directory}{agn_filename}',dpi = 300)
    plt.show()
    plt.close()
    '''
    ## plotting all of the morph_clumps in array
    # plot the fraction of disks that have clumps to histogram
    c=0
    for i in morph_clumps:
        plt.ylabel('Count')
        plt.xlabel('AGN Fraction (%)')
        plt.hist(morph_agn[c], label=labels[c],bins=bins,lw=1,fill=False,histtype='step',color=bin_colors[labels[c]])
        plt.hist(i['AGN(%)'],bins=bins,fill=True,label='Has Clumps')
        plt.title(f"AGN Fraction by Galaxy Type ({labels[c]} w/ Clumps)")
        plt.legend()
        plt.savefig(f'{working_directory}agn-hist-w-clumps-{labels[c]}',dpi = 300)
        # plt.show()
        plt.close()
        c+=1
    ## plotting all of the morph_spiral in array
    # plot the fraction of disks that have spiral arms to histogram
    c=0
    for i in morph_spiral:
        plt.ylabel('Count')
        plt.xlabel('AGN Fraction (%)')
        plt.hist(morph_agn[c], label=labels[c],bins=bins,lw=1,fill=False,histtype='step',color=bin_colors[labels[c]])
        plt.hist(i['AGN(%)'],bins=bins,fill=True,label='Has Arms')
        plt.title(f"AGN Fraction by Galaxy Type ({labels[c]} w/ Spiral Arms)")
        plt.legend()
        plt.savefig(f'{working_directory}agn-hist-w-spirals-{labels[c]}',dpi = 300)
        # plt.show()
        plt.close()
        c+=1
    '''

    # ### get the names and display cutouts of the jades rgb images for different ranges of each
    # irr = irregular_only.sort_values(by='AGN(%)')
    # print('IRREGULARS -----')
    # print(irr[['ID','AGN(%)']]) # get both properties & agn frac
    # sph = spheroid_only.sort_values(by='AGN(%)')
    # print('SPHEROIDS -----')
    # print(sph[['ID','AGN(%)']]) # get both properties & agn frac
    # disk = disk_only.sort_values(by='AGN(%)')
    # print('DISKS -----')
    # print(disk[['ID','AGN(%)']]) # get both properties & agn frac


def type_z_hist(): #### need to update this one if we're gonna keep using it ** (7/15)
    #### make a hist to showcase the type of galaxy and redshift distribution of the data
    histbins = 5
    plt.title("Data Distribution of Redshifts (z)")
    plt.ylabel('Count')
    plt.xlabel('z')
    hist_bins = np.histogram(np.hstack((morph_z[0]['z'],morph_z[1]['z'],morph_z[2]['z'],morph_z[3]['z'])),bins=histbins)[1] # ideally would put all them but I'm lazy and it looks good as is
    for i in morph_z:
        plt.hist(i['z'],bins=hist_bins,lw=1,fill=False,histtype='step')
    plt.legend(labels,prop={"size":8})
    plt.savefig(f'{working_directory}z_distribution_type',dpi=300)
    plt.show()

    # bins=np.histogram(np.hstack((disk_z,sph_z,irr_z,disk_sph_z,disk_irr_z,disk_irr_sph_z,orph_z)), bins=histbins)[1] #get the bin edges


def agn_frac_z(): ##### need to update this one
    #### make a scatter plot of the agn fraction vs redshift to showcase if there is a relation there?


    ### Calculate the % morphology in each agn section (<20%, 20<x<80%, >80%)

    ### old but save for later bc this shorthand is super OP
    # sf_sample = [i.query('`AGN(%)`<=20') for i in morph_types]
    # comp_sample = [i.query('`AGN(%)`>20 & `AGN(%)`<80') for i in morph_types]
    # pure_agn_sample = [i.query('`AGN(%)`>=80') for i in morph_types]
    # Create GridSpec layout
    fig = plt.figure(figsize=(14, 6))
    gs = GridSpec(3, 2, figure=fig)
    ax_scatter = fig.add_subplot(gs[:, 0])
    ax_pie_1 = fig.add_subplot(gs[0, 1])
    ax_pie_2 = fig.add_subplot(gs[1, 1], position=[0.56, 0.55, 0.18, 0.35])  # Adjust position for better fit
    ax_pie_3 = fig.add_subplot(gs[2, 1], position=[0.75, 0.55, 0.18, 0.35])  # Adjust position for better fit

    # ax_clumps = fig.add_subplot(gs[0,2])

    # fig, ax = plt.subplots()
    # fig.tight_layout(rect=[0,0,.8,.9])
    # ftsz = 8




    sf_sample = np.empty(len(morph_types),dtype=object)
    comp_sample = np.empty(len(morph_types),dtype=object)
    pure_agn_sample = np.empty(len(morph_types),dtype=object)
    sf_nums,comp_nums,pure_agn_nums = [np.zeros(len(morph_types)) for i in range(3)]

    ### 3 lines below replaced with table.query
    sf_sample_table = combined.query('`AGN(%)`<20')
    comp_sample_table = combined.query('`AGN(%)`>=20 & `AGN(%)`<80')
    pure_sample_table = combined.query('`AGN(%)`>=80')
    total_sf = sf_sample_table.shape[0]
    total_comp = comp_sample_table.shape[0]
    total_pure_agn = pure_sample_table.shape[0]

    v=0
    for i in morph_types:
        # plotting the scatter
        bin_type = list(bin_colors)[v]
        color = list(bin_colors.values())[v]
        ax_scatter.scatter(i['z'],i['AGN(%)'],color=color, label=bin_type)

        # getting the elements
        sf_sample[v] = i.query('`AGN(%)`<=20')
        comp_sample[v] = i.query('`AGN(%)`>20 & `AGN(%)`<80')
        pure_agn_sample[v] = i.query('`AGN(%)`>=80')



        sf_nums[v] = (sf_sample[v].shape[0]/total_sf)*100

        comp_nums[v] = (comp_sample[v].shape[0]/total_comp)*100
        pure_agn_nums[v] = (pure_agn_sample[v].shape[0]/total_pure_agn)*100
        v+=1


    # Pie chart for AGN fraction >= 80%
    colors = bin_colors.values()
    ax_pie_1.pie(pure_agn_nums, autopct='%1.0f%%', startangle=140, colors=colors)
    has_clumps_percent = (pure_sample_table.query('`Has Clumps (flag)`=="yes"').shape[0]/total_pure_agn)*100
    has_spiral_arms_percent = (pure_sample_table.query('`Spiral Arms (flag)`=="yes"').shape[0]/total_pure_agn)*100
    is_merger_percent = (pure_sample_table.query('`Merger (flag threshold)`=="yes"').shape[0]/total_pure_agn)*100
    ax_pie_1.text(1.2, 0.5, f'Has Clumps: {has_clumps_percent:.0f}%\nHas Spiral Arms: {has_spiral_arms_percent:.0f}%\nMerger: {is_merger_percent:.0f}%', transform=ax_pie_1.transAxes)
    ax_pie_1.set_title('Pure AGN Distributions')

    # Pie chart for AGN fraction between 20 & 80%
    ax_pie_2.pie(comp_nums, autopct='%1.0f%%', startangle=140, colors=colors)
    has_clumps_percent = (comp_sample_table.query('`Has Clumps (flag)`=="yes"').shape[0]/total_comp)*100
    has_spiral_arms_percent = (comp_sample_table.query('`Spiral Arms (flag)`=="yes"').shape[0]/total_comp)*100
    is_merger_percent = (comp_sample_table.query('`Merger (flag threshold)`=="yes"').shape[0]/total_comp)*100
    ax_pie_2.text(1.2, 0.5, f'Has Clumps: {has_clumps_percent:.0f}%\nHas Spiral Arms: {has_spiral_arms_percent:.0f}%\nMerger: {is_merger_percent:.0f}%', transform=ax_pie_2.transAxes)
    ax_pie_2.set_title('Composite Distributions')

    # Pie chart for AGN fraction <= 20%
    ax_pie_3.pie(sf_nums, autopct='%1.0f%%', startangle=140, colors=colors)
    has_clumps_percent = (sf_sample_table.query('`Has Clumps (flag)`=="yes"').shape[0]/total_sf)*100
    has_spiral_arms_percent = (sf_sample_table.query('`Spiral Arms (flag)`=="yes"').shape[0]/total_sf)*100
    is_merger_percent = (sf_sample_table.query('`Merger (flag threshold)`=="yes"').shape[0]/total_sf)*100
    ax_pie_3.text(1.2, 0.5, f'Has Clumps: {has_clumps_percent:.0f}%\nHas Spiral Arms: {has_spiral_arms_percent:.0f}%\nMerger: {is_merger_percent:.0f}%', transform=ax_pie_3.transAxes)
    ax_pie_3.set_title('SF Distributions')

    plt.tight_layout(pad=2.)
    # # Adding labels and title for the scatter plot
    # ax_scatter.set_xlabel('Redshift')
    # ax_scatter.set_ylabel('AGN Fraction')
    # ax_scatter.set_title('Galaxy Data Scatter Plot')
    # ax_scatter.legend(loc='best')
    # plt.show()
    # exit()
    '''
    v=0
    sf_len=comp_len=pure_len=clump_len=spiral_len=merger_len=0
    for i in morph_types:
        #### get the fractional counts from the classifications and flags for each sub section (sf, composite source, pure agn) ####

        #### for the main classifications (disk, spheroid, irregular, combinations of them...)
        sf_num = (sf_sample[v].shape[0]/total_sf)*100
        comp_num = (comp_sample[v].shape[0]/total_comp)*100
        pure_agn_num = (pure_agn_sample[v].shape[0]/total_pure_agn)*100
        ax.set(xlim=[0, 5])
        if(sf_num>0):
            print(f'SF % for {labels[v]}: {sf_num:.2f} %')
            ax.annotate(f'{labels[v]}:{sf_num:.0f}%',xy=(3.5,sf_len),xytext=(4,19-sf_len),horizontalalignment='center',verticalalignment='top',clip_on=False,fontsize=ftsz)
            sf_len+=4
        if(comp_num>0):
            print(f'Composite % for {labels[v]}: {comp_num:.2f} %')
            ax.annotate(f'{labels[v]}:{comp_num:.0f}%',xy=(3.5,comp_len),xytext=(4,70-comp_len),horizontalalignment='center',verticalalignment='top',clip_on=False,fontsize=ftsz)
            comp_len+=4
        if(pure_agn_num>0):
            print(f'Pure AGN % for {labels[v]}: {pure_agn_num:.2f} %')
            ax.annotate(f'{labels[v]}:{pure_agn_num:.0f}%',xy=(3.5,pure_len),xytext=(4,98-pure_len),horizontalalignment='center',verticalalignment='top',clip_on=False,fontsize=ftsz)
            pure_len+=4

        print()
        v+=1
    '''

    '''
    #### now for the flags (clumps, spiral arms, and mergers)
    # sf_clumps = [f.query('`Has Clumps (flag)`=="yes"') for f in sf_sample]
    # sf_clump_num = np.sum([f.shape[0] for f in sf_clumps]) # what this is awesome! so easy, so inefficient!
    # comp_clumps = [f.query('`Has Clumps (flag)`=="yes"') for f in comp_sample]
    # comp_clump_num = np.sum([f.shape[0] for f in comp_clumps])
    # pure_clumps = [f.query('`Has Clumps (flag)`=="yes"') for f in pure_agn_sample]
    # pure_clump_num = np.sum([f.shape[0] for f in pure_clumps])

    # sf_clump = sf_clump_num/total_sf*100
    # comp_clump = comp_clump_num/total_comp*100
    # pure_clump = pure_clump_num/total_pure_agn*100

    # if(sf_clump>0): ax.annotate(f'Has Clumps: {sf_clump:.0f}%',xy=(5.25,0),xytext=(5.2,18),horizontalalignment='left',verticalalignment='top',annotation_clip=False,fontsize=ftsz)
    # if(comp_clump>0): ax.annotate(f'Has Clumps: {comp_clump:.0f}%',xy=(5.25,0),xytext=(5.2,70),horizontalalignment='left',verticalalignment='top',annotation_clip=False,fontsize=ftsz)
    # if(pure_clump>0): ax.annotate(f'Has Clumps: {pure_clump:.0f}%',xy=(5.25,0),xytext=(5.2,95),horizontalalignment='left',verticalalignment='top',annotation_clip=False,fontsize=ftsz)

    # clumps
    samples = np.array([sf_sample,comp_sample,pure_agn_sample],dtype=object)
    totals = np.array([total_sf,total_comp,total_pure_agn])
    ys = np.array([18,70,95])
    a=0
    for s in samples:
        clump_total = np.sum([i.query('`Has Clumps (flag)`=="yes"').shape[0] for i in s])
        clumps=clump_total/totals[a]*100
        if(clumps>0): ax.annotate(f'Has Clumps: {clumps:.0f}%',xy=(5.25,0),xytext=(5.2,ys[a]),horizontalalignment='left',verticalalignment='top',annotation_clip=False,fontsize=ftsz)
        a+=1


    # spirals
    ys = np.array([13,65,90])
    a=0
    for s in samples:
        spiral_total = np.sum([i.query('`Spiral Arms (flag)`=="yes"').shape[0] for i in s])
        spiral=spiral_total/totals[a]*100
        if(spiral>0): ax.annotate(f'Has Spiral Arms: {spiral:.0f}%',xy=(5.25,0),xytext=(5.2,ys[a]),horizontalalignment='left',verticalalignment='top',annotation_clip=False,fontsize=ftsz)
        a+=1

    # merger fraction
    ys = np.array([8,60,85])
    a=0
    for s in samples:
        merger_total = np.sum([i.query('`Merger (flag threshold)`=="yes"').shape[0] for i in s])
        merger=merger_total/totals[a]*100
        if(merger>0): ax.annotate(f'Mergers: {merger:.0f}%',xy=(5.25,0),xytext=(5.2,ys[a]),horizontalalignment='left',verticalalignment='top',annotation_clip=False,fontsize=ftsz)
        a+=1



    '''

    # plt.legend(labels,loc='best',prop={'size': 7})
    ax_scatter.legend(labels,loc='upper left',prop={'size':8})
    ax_scatter.set_xlabel('z')
    ax_scatter.set_ylabel('AGN(%)')
    ax_scatter.axhline(20,linestyle='--',color='lightgray')
    ax_scatter.axhline(80,linestyle='--',color='lightgray')
    ax_scatter.set_title('AGN(%) & Redshift(z) of Galaxy Types')

    plt.savefig(f'{working_directory}agn_fraction_z_pie_chart',dpi=500)
    plt.show()



def agn_frac_merger_hist(): ### update this is we still want to use it*** (7/15)
    #### make a histogram of the agn fraction of mergers (then add cutouts later)
    #### to show if there's a relation bewteen the visual classification and quantitative measurements (Gini/M20) (to be found later)
    histbins = 10
    hist_bins = np.histogram(np.hstack((morph_types[0]['AGN(%)'],morph_types[1]['AGN(%)'],morph_types[2]['AGN(%)'],morph_types[3]['AGN(%)'])),bins=histbins)[1]
    # plt.hist(mergers['AGN(%)'],bins=hist_bins,lw=1.5,fill=False,color=bin_colors[mergers['Classification']],histtype='step')
    merger_data = mergers.groupby(['AGN(%)','Classification']).size().unstack(fill_value=0)
    merger_data.plot(kind='bar',stacked=True, figsize=(10,6))
    # Highlight the agn regions
    y_max = plt.gca().get_ylim()[1]
    x_max = plt.gca().get_xlim()[1]

    # plt.axvspan(80, 100, color='lightgray', alpha=0.2)
    # plt.text(x=.90*x_max, y=.9*y_max, s='AGN', fontsize=12, ha='center', va='center')
    # plt.axvspan(20, 80, color='lightgray', alpha=0.2)
    # plt.text(x=.50*x_max, y=.9*y_max, s='Composite', fontsize=12, ha='center', va='center')
    # plt.axvspan(-5, 20, color='lightgray', alpha=0.2)
    # plt.text(x=.1*x_max, y=0.9*y_max, s='SF', fontsize=12, ha='center', va='center')
    # mergers['AGN(%)'].hist()
    # for i in mergers:
    #     print(i)
    #     plt.hist(i['AGN(%)'],bins=hist_bins,lw=1.5,fill=False,color=bin_colors[i['Classification'][0]],histtype='step')
    plt.legend(labels,loc='best',prop={'size': 8})
    plt.xlabel('AGN(%)')
    plt.ylabel('Count')

    plt.title('AGN(%) of Mergers')
    plt.savefig(f'{working_directory}agn_merger_hist',dpi=300)
    plt.show()
    plt.close()

def agn_frac_merger_scatter(): ### update this is we still want to use it*** (7/15)
    #### make a scatter plot of the agn fraction vs merger (% agreement) with color coded classifications
    # for i in mergers:
        # plt.scatter(i['AGN(%)'],i['Merger? (%)'])
    merger_data = mergers.groupby(['AGN(%)','Classification']).size().unstack(fill_value=0)
    merger_data.plot(kind='scatter',figsize=(10,6))
    plt.scatter(mergers['AGN(%)'],mergers['Merger? (%)'])
    plt.xlabel('AGN(%)')
    plt.ylabel('Agrement on Merger (%)')
    # plt.xlim([0,50])
    # plt.ylim([0,100])
    plt.axhline(61,linestyle='--',color='darkgray',label='Yes Threshold')
    # plt.axhline(66/2,linestyle='--',color='lightgray')
    # labels.append('Yes Threshold')
    # labels.append('Maybe Threshold')
    plt.title('AGN(%) of Mergers')
    # plt.legend(labels,loc='upper right',prop={'size': 7})
    plt.legend()
    plt.savefig(f'{working_directory}agn_merger_scatter',dpi=300)
    plt.show()
    plt.close()
    # labels.remove('Yes Threshold')
    # labels.remove('Maybe Threshold')

def compare_subsample():
    # compare the sub sample classifications of the data
    # get the last time the file was updated
    ti_m = os.path.getmtime(detailed_file)
    m_ti = time.ctime(ti_m)
    print(f"[!!!] Detailed file ({details_filename}) was last updated at: {m_ti}")

    ### import the data from the file name near the imports
    df = pd.read_csv(detailed_file,sep='\t')
    df.sort_values(by='ID',ascending=True,inplace=True)
    # remove_excluded(df,name='ID')
    # sub_sil = df.query('`AGN^a`>=80')
    # sub_comp = df.query('`AGN^a`<80 & `AGN^a`>20')
    
    # table.insert
    combined.insert(2,"Sub-sample",df['Sub-sample'].values)
    # table.groupby
    count_data = combined.groupby(['Classification', 'Sub-sample']).size().unstack(fill_value=0)
    count_data.plot(kind='barh',stacked=True, figsize=(10,6))
    plt.tight_layout(pad=2.)
    plt.title('Sub-sample Counts by Classification Type')
    plt.savefig(f'{working_directory}subsample_counts',dpi=300)
    plt.show()
    plt.close()

    '''
    sub_sil = df.query('`Sub-sample`=="Sil AGN"')[['ID','AGN^a']]
    sub_feat = df.query('`Sub-sample`=="Feat AGN"')[['ID','AGN^a']]
    sub_unknown = df.query('`Sub-sample`=="..."')[['ID','AGN^a']]

    # frame = table.where(table['ID']==sub_sil)
    # frame = table[(table.Properties == sub_sil)]
    sub_sil_classify = table[(table['ID'].isin(sub_sil['ID']))].sort_values(by='ID',ascending=True)
    sub_feat_classify = table[(table['ID'].isin(sub_feat['ID']))].sort_values(by='ID',ascending=True)
    sub_unknown_classify = table[(table['ID'].isin(sub_feat['ID']))].sort_values(by='ID',ascending=True)

    ### here we need to get the individual classifications from the filter method, then plot a histogram or something
    ### across the agn fraction on x axis
    sil_class,sil_labels = add_classifications(sub_sil_classify)
    comp_class,comp_labels = add_classifications(sub_feat_classify)
    c=0
    for i in sil_class:
        plt.hist(i['AGN(%)'],bins=5)
        # plt.hist(comp_class[c]['AGN(%)'],bins=5)
        c+=1
    # for i in comp_class:
        # plt.hist(i['AGN(%)'])
    plt.legend(labels=comp_labels)
    # [plt.hist(i['AGN(%)'] for i in comp_class)]
    # plt.hist(comp_class[0]['AGN(%)'],bins=5)
    # plt.hist(sub_comp_classify['AGN(%)'],bins=5)
    plt.show()
    # # plt.axhline(80,linestyle='--',color='lightgray',label='Sil AGN Threshold')
    # plt.axvspan(0,3,ymin=.80,ymax=1,alpha=0.15,color='gray',label='Sil AGN')
    # plt.legend(loc='best')
    # plt.xlabel('z')
    # plt.ylabel('AGN(%)')
    # plt.title('AGN(%) & Z of Subsamples')
    # plt.show()
    '''

### done with the visual comparison, moving to statmorph analysis


#### already have the statmorph files imported *******



### testing the problem sources
#### searching for point sources
# point_sources = ['GS_IRS25','GN_IRS55']
# xsearch = 'M20'
# ysearch = 'Gini'

# for id in point_sources:
#     source_row = combined.loc[combined['ID'] == id]
#     plt.scatter(source_row[xsearch],source_row[ysearch])
#     plt.annotate(source_row['ID'].iloc[0],(source_row[xsearch],source_row[ysearch]))
    
#  ## show the Lotz et al. 2008 merger line
#     x = [0,-3.0]
#     y = [.33,.75]
#     plt.plot(x, y,color='darkgray',label='Lotz+2008')
#     plt.annotate('Merger',(-.75,0.65),color='darkgray')

#     ## show E/S0/Sa & Sb/Sc/Ir separator line (Lotz et al. 2008)
#     x = [-1.7,-3.0]
#     y = [0.568,0.38]
#     plt.plot(x, y,color='darkgray')
#     plt.annotate('E/S0/Sa',(-2.5,0.55),color='darkgray')
#     plt.annotate('Sb/Sc/Ir',(-1.2,0.42),color='darkgray')

# plt.xlim([1.8,4.6])
# plt.ylim([0,0.85])
# plt.gca().invert_yaxis()
# plt.xlabel(xsearch)
# plt.ylabel(ysearch)
# plt.show()
# plt.close()
# exit()


#### testing with the KS test to see if anything alerting?
# from scipy.stats import kstest
# data = combined[combined['AGN(%)']>=80]['Asymmetry (A)']
# # data = combined[combined['AGN(%)']>=80]['Concentration (C)']
# # data = combined[combined['AGN(%)']>=80]['Gini']
# # data = combined[combined['AGN(%)']>=80]['M20']

# # Perform the K-S test against the normal distribution
# statistic, p_value = kstest(data, 'norm')

# # Print the results
# print(f"KS statistic: {statistic}")
# print(f"P-value: {p_value}")
# exit()
###########################################################


## showing k-correction evolution of single source
def k_corr_evolution():
    # read values from file
    # df = pd.read_csv('/home/robbler/research/making_psf_edits/testing_k_corr.txt', delimiter=':')
        
    file = open('/home/robbler/research/making_psf_edits/testing_k_corr.txt','r')
    data = file.read()
    file.close()

    # Parse the data
    lines = data.strip().split('\n')
    parsed_data = []
    temp_dict = {}
    for line in lines:
        if line.startswith('Source:'):
            if temp_dict:
                parsed_data.append(temp_dict)
                temp_dict = {}
            temp_dict['Source'] = line.split(':')[1]
        elif line.startswith('Filter:'):
            temp_dict['Filter'] = line.split(':')[1].strip()
        elif line.startswith('Asymmetry:'):
            temp_dict['Asymmetry'] = float(line.split(':')[1].strip())
        elif line.startswith('Gini:'):
            temp_dict['Gini'] = float(line.split(':')[1].strip())

    # Append the last entry
    if temp_dict:
        parsed_data.append(temp_dict)

    # Create a DataFrame
    df = pd.DataFrame(parsed_data)

    # Display the DataFrame
    print(df)

    # Merge the data on the Source column
    merged_df = df.groupby('Source').apply(lambda x: x[['Filter', 'Asymmetry', 'Gini']].to_dict(orient='list')).reset_index()
    merged_df.columns = ['Source', 'Details']


    agn_obs = {
        'GN_IRS11':['23',2.776],
        'GN_IRS17':['47',3.556],
        'GN_IRS27':	['26',1.988],
        'GS_IRS14':	['50',2.776],
        'GS_IRS25':	['0',1.988],
        'GS_IRS60':['40',2.776],
        'GN_IRS2':	['0',1.988],
        'GS_IRS62':	['0',1.501],
        'GS_IRS9': ['80(psf-effects)',4.401]
           }
    # Display the merged DataFrame
    # print(merged_df[['Source','Details']])
    fig,ax = plt.subplots(2,1)
    filter_options = [1.501,1.988,2.776,3.566,4.401]
    for id in merged_df.iterrows():
        source = id[1]['Source'].strip()
        obs_filt = agn_obs[source][1]
        agn = agn_obs[source][0]
        xs = filter_options
        y_a = id[1]['Details']['Asymmetry']
        y_g = id[1]['Details']['Gini']
        fig,ax = plt.subplots(2,1)
        ax[0].plot(xs,y_a,'o-')
        ax[0].set_ylabel('Asymmetry')
        ax[0].axvline(obs_filt,linestyle='--')
        ax[1].axvline(obs_filt,linestyle='--')
        ax[1].plot(xs,y_g,'o-')
        ax[1].set_ylabel('Gini')
        plt.xlabel('Wavelength ($\mu m$)')
        plt.suptitle(f'{source} (.{agn} agn) K-correction Idea')
        plt.savefig(f'/home/robbler/research/making_psf_edits/{source}_k_corr_a_gini.png',dpi=200)
        # plt.show()

def conv_stat_vals():
    sources = ['GN_IRS36(.93)','GN_IRS43(.84)','GS_IRS9(.80)','GN_IRS4 (.17)']
    c_s = [[1.772148129,1.583310724,4.341623673,2.698351791],[1.772148129,1.583310724,4.341623673,2.698351791],[4.324636132,2.854352176,3.064040347,2.675849208],[4.324636132,2.854352176,3.064040347,2.675849208]]
    a_s = [[0.075922317,0.937637575,0.064693479,0.211675086],[0.075922317,0.937637575,0.064693479,0.211675086],[0.151905164,0.066915952,0.040499877,0.275871043],[0.151905164,0.066915952,0.040499877,0.275871043]]
    ginis = [[0.312732191,0.460723916,0.324699901,0.511568396],[0.312732191,0.460723916,0.324699901,0.511568396],[0.601423313,0.584647868,0.530283013,0.513032632],[0.601423313,0.584647868,0.530283013,0.513032632]]
    m20s = [[-1.573060311,-1.469682853,-1.366981295,-1.695406395],[-1.573060311,-1.469682853,-1.366981295,-1.695406395],[-2.274841868,-1.702513933,-1.786607639,-1.637080716],[-2.274841868,-1.702513933,-1.786607639,-1.637080716]]

    x_key = ['Conv & PSF','Conv no PSF','PSF no Conv','No PSF or Conv']
    fig,ax = plt.subplots(2,2)
    for i in range(len(sources)):
        concs = [c_s[f][i] for f in range(4)]
        ax[0,0].plot(x_key,concs,linestyle='-',marker='o',label=f'{sources[i]}')
        ax[0,0].set_ylabel('Concentration')

        asym = [a_s[f][i] for f in range(4)]
        ax[1,0].plot(x_key,asym,linestyle='-',marker='o',label=f'{sources[i]}')
        ax[1,0].set_ylabel('Asymmetry')

        gini_vals = [ginis[f][i] for f in range(4)]
        ax[1,1].plot(x_key,gini_vals,linestyle='-',marker='o',label=f'{sources[i]}')
        ax[1,1].set_ylabel('Gini')

        m20_vals = [m20s[f][i] for f in range(4)]
        ax[0,1].plot(x_key,m20_vals,linestyle='-',marker='o',label=f'{sources[i]}')
        ax[0,1].set_ylabel('M20')
        # ax[0,0].scatter(c_s[0][i],a_s[0][i],label=f'{sources[i]}')
        # ax[0,0].set_title('Convolution & PSF inputted')
        # ax[0,0].set_xlabel('Concentration')
        # ax[0,0].set_ylabel('Asymmetry')
        # ax[1,0].set_xlabel('Concentration')
        # ax[1,0].set_ylabel('Asymmetry')
        # ax[1,1].set_xlabel('Concentration')
        # ax[1,1].set_ylabel('Asymmetry')
        # ax[0,1].set_xlabel('Concentration')
        # ax[0,1].set_ylabel('Asymmetry')

        # ax[1,0].scatter(c_s[1][i],a_s[1][i],label=f'{sources[i]}')
        # ax[1,0].set_title('Convolution & without PSF input')

        # ax[0,0].scatter(c_s[2][i],a_s[2][i],label=f'{sources[i]}')
        # ax[0,0].set_xlabel('Convolution & PSF inputted')
    plt.legend()
    plt.show()

# conv_stat_vals()
# exit()

## S/N histogram
def s_n_hist():
    bin_edges = [0,10,15,20,30,40,50,70,100,150,200,250] # Smaller bins for smaller values
    plt.hist(combined['S/N'],bins=bin_edges,edgecolor='black')
    med = np.median(combined['S/N'])
    std = np.std(combined['S/N'])
    plt.axvline(med,linestyle='--',c='red',label=f'$\\tilde{{S/N}}={med:.1f}$')
    plt.xlim(0,50)
    plt.xticks([0,5,10,15,20,25,30,35,40,45,50])
    # plt.axvline(med+std,linestyle='--',c='red')
    # plt.axvline(med-std,linestyle='--',c='red')
    plt.title('S/N of Sources')
    plt.xlabel('S/N')
    plt.legend()
    plt.show()

# s_n_hist()
# exit()

##### Error bar plots #####
def test_plot_error_bars():
    ### for GN_IRS10
    ### Gaussian noise applied scaled by 1/SN = 1/~50
    conc = np.array([3.0550657057031216,0.002794759982496231, 0.0030080543239394686])
    asym = np.array([0.2621439748005437,0.0033458888329225, 0.0018836116691538785])
    gini = np.array([0.5067456027053324,0.0016034678642217193, 0.0019597491132878897])
    m20 = np.array([-1.791144591569481,0.00299742678402537, 0.0034045452864090997])
    sn = 53.55214309692383
    agn_perc = 0


    ### for 3 lowest S/N galaxies
    # conc_1 = np.array([3.866892754, 0.022605117, 0.027190218])
    # conc_2 = np.array([3.299109852,0.010404236,0.008720217])
    # conc_3 = np.array([3.660702562,0.014971688,0.015575261])

    conc_1 = np.array([3.660702562,0.005830367,0.004887749])
    conc_2 = np.array([3.045572267,0.001298482,0.001208933])
    conc_3 = np.array([3.552820596,0.00070891,0.000910829])
    agn_perc = np.array([0,100,66])
    sn = np.array([18.85951424,20.39265442,17.11173248])

    # plt.plot(agn_perc, conc[0],'bo')
    plt.title(f'3 lowest sn galaxies (SN AVG: 18)')
    plt.ylim([3.0,4.0])
    # plt.xlim([-1,100])
    c=0
    for i in enumerate([conc_1,conc_2,conc_3]):
        plt.errorbar(agn_perc[i[0]],i[1][0],yerr=i[1][1:2],fmt='.')
        c+=1
    # print(f'SN:{sn:.2f}')
    plt.show()
    exit()

# test_plot_error_bars()

def error_sn_plot():
    ### plot the relative error for each respective measurement
    # with respect to the S/N to see it decrease with increasing S/N    

    ### plot all 4 measurements we're focued on and have error for (Concentration, asymmetry, gini, m20)
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))

    filtered_group = combined.groupby('Filter')
    for i in filtered_group:
        # concentration first
        ax[0,0].plot(i[1]['S/N'],(i[1]['Concentration (C) Error (16%)']+i[1]['Concentration (C) Error (84%)'])/2,'.',label=f'{i[0]}')
        # asymmetry
        ax[0,1].plot(i[1]['S/N'],(i[1]['Asymmetry (A) Error (16%)']+i[1]['Asymmetry (A) Error (84%)'])/2,'.',label=f'{i[0]}')
        # gini
        ax[1,0].plot(i[1]['S/N'],(i[1]['Gini Error (16%)']+i[1]['Gini Error (84%)'])/2,'.',label=f'{i[0]}')
        # m20
        ax[1,1].plot(i[1]['S/N'],(i[1]['M20 Error (16%)']+i[1]['M20 Error (84%)'])/2,'.',label=f'{i[0]}')



    # ax[0,0] = concentration
    ax[0,0].legend()
    # ax[0,0].set_title('Concentration error (*100) vs S/N by Filter')
    ax[0,0].set_ylabel('Concentration (C) Error')
    ax[0,0].set_yscale('log')
    ax[0,0].set_xlabel('S/N')
    ax[0,0].set_xscale('log')
    # plt.xlim([0,100])

    # ax[0,1] = asymmetry
    # ax[0,1].legend()
    # ax[0,0].title('Asymmetry error (*100) vs S/N by Filter')
    ax[0,1].set_ylabel('Asymmetry (A) Error')
    ax[0,1].set_yscale('log')
    ax[0,1].set_xscale('log')
    # ax[0,1].set_xlabel('S/N (avg per pixel)')

    # ax[1,0] = gini
    # ax[1,0].legend()
    ax[1,0].set_ylabel('Gini Error')
    ax[1,0].set_yscale('log')
    ax[1,0].set_xscale('log')
    # ax[1,0].set_xlabel('S/N (avg per pixel)')

    # ax[0,1] = asymmetry
    # ax[1,1].legend()
    ax[1,1].set_ylabel('M20 Error')
    # ax[1,1].set_yscale('log')
    # ax[1,1].set_xscale('log')
    # ax[1,1].set_xlabel('S/N (avg per pixel)')

    # ax[0,0].set_xlim([0,100])
    # ax[1,0].set_xlim([0,100])
    # ax[1,1].set_xlim([0,100])
    # ax[0,1].set_xlim([0,100])
    # plt.tight_layout()
    fig.tight_layout(pad=2)
    fig.suptitle('Avg Relative Error vs S/N (per Filter)',fontsize=16)
    plt.show()

# error_sn_plot()
# exit()

def error_dist_hist():
    ### relative error distribution for each filter in histogram form
    ### we want to see a gaussian distribution here for each measurement
    ### not focusing on S/N, just what the overall distribution of the values are
    ### this should be a gaussian distribution if statmorph isn't messing with it
    ### if it isn't gaussian we need to adjust the error bars to be different than the 16th & 84th percentile & median for root value
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))

    filtered_group = combined.groupby('Filter')

    for i in filtered_group:
        ### for now it's just the histogram of a few sources, will eventually do the full thing but it'll take a while to run all them and save them
        ### idk if I can save an np array to a table like that        
        if(i[0]=='f150w'):
            bins=1 
        if(i[0]=='f200w'):
            bins=1
        if(i[0]=='f277w'):
            bins=8
        if(i[0]=='f356w'):
            bins=1
        if(i[0]=='f444w'):
            bins=1 
        # ax[0,0].hist(conc-np.median(conc),5,alpha=0.5,label=i[0])
        ax[0,0].hist((i[1]['Concentration (C) Error (84%)']+i[1]['Concentration (C) Error (16%)'])/2,bins=bins,label=i[0],histtype='step',linewidth=2)
        # # # plt.xlim([-10,10])
        ax[0,1].hist((i[1]['Asymmetry (A) Error (84%)']+i[1]['Asymmetry (A) Error (16%)'])/2,bins=bins,label=i[0],histtype='step',linewidth=2)
        # ax[1,0].hist(asym-np.median(asym),5,alpha=0.5,label=i[0])
        ax[1,0].hist((i[1]['Gini Error (84%)']+i[1]['Gini Error (16%)'])/2,bins=bins,label=i[0],histtype='step',linewidth=2)
        # ax[1,1].hist(gini-np.median(gini),5,alpha=0.5,label=i[0])
        ax[1,1].hist((i[1]['M20 Error (84%)']+i[1]['M20 Error (16%)'])/2,bins=2*bins,label=i[0],histtype='step',linewidth=2)
        # ax[0,1].hist(m20-np.median(m20),5,alpha=0.5,label=i[0])


    ax[0,0].set_xlabel('Concentration')
    ax[0,0].legend()
    ax[0,1].set_xlabel('Asymmetry')
    ax[1,0].set_xlabel('Gini')
    ax[1,1].set_xlabel('M20')
    fig.suptitle('Avg Error Histogram by Filter')
    ax[0,0].set_xscale('log')
    ax[0,1].set_xscale('log')
    ax[1,0].set_xscale('log')
    ax[1,1].set_xscale('log')
    plt.show()

# error_dist_hist()
# exit()

def error_hist_for_source(id='GS_IRS12'):
    ### for plotting the distribution of values for each, just for testing purposes
        fig, axes = plt.subplots(2, 2, figsize=(10, 10))

        source_row = combined.loc[combined['ID'] == id]
        c_err = import_data.read_errorbar_data(source_row['Concentration (C) Error Full'])
        a_err = import_data.read_errorbar_data(source_row['Asymmetry (A) Error Full'])
        g_err = import_data.read_errorbar_data(source_row['Gini Error Full'])
        m_err = import_data.read_errorbar_data(source_row['M20 Error Full'])
        labels = ['Concentration', 'Asymmetry', 'Gini', 'M20']
        data = [c_err, a_err, g_err, m_err]
        ### for GS_IRS12 plot the root values with no noise
        if(id=='GS_IRS12'):
            root_vals = [3.0455722671575853,0.17032136048318505,0.5772125623193657,-1.4582065574879546]
        elif(id=='GN_IRS17'):
            # morph.concentration=2.7793204937496485
            # morph.asymmetry=0.234235029703258
            # morph.gini=0.4785028770680358
            # morph.m20=-1.612607324410134
            root_vals = [2.7793204937496485,0.234235029703258,0.4785028770680358,-1.612607324410134]

        axes = axes.flatten()
        for i in range(4):
            axes[i].hist(data[i], bins=20, color='b', alpha=0.7)
            axes[i].set_title(labels[i])
            axes[i].axvline(np.median(data[i]), color='r', linestyle='dashed', linewidth=1.5)
            axes[i].axvline(root_vals[i],color='g',linestyle='dashed',linewidth=1.5)
            fig.suptitle(f"Distribution of Values for {id}", fontsize=16)

        # savepath = (f'research/statmorph_output/grizli/{search_id}_Error_Histogram')
        # plt.savefig(savepath,dpi=100)
        plt.show()
        plt.close()

# error_hist_for_source('GN_IRS17')
# exit()

def sn_on_error(type='combined'):
    ### three types: 
    # ---- separate (in four separate histograms)
    # ---- combined (all in one histogram)
    # function for the s/n on the error of each measurement
    # should give us an idea of how much we could trust our measurements


    bins = 20
    # print(i[1]['ID'])
    # exit()
    # ax[0,0].hist(conc-np.median(conc),5,alpha=0.5,label=i[1])
    if(type=='separate'):
        fig, ax = plt.subplots(2,2, figsize=(10, 8))
        ax[0,0].hist(combined['Concentration (C)']/((combined['Concentration (C) Error (84%)']+combined['Concentration (C) Error (16%)'])/2),label='Concentration',bins=bins,histtype='step',linewidth=2)
        ax[0,1].hist(combined['Asymmetry (A)']/((combined['Asymmetry (A) Error (84%)']+combined['Asymmetry (A) Error (16%)'])/2),label='Asymmetry',bins=bins,histtype='step',linewidth=2)
        ax[1,0].hist(combined['Gini']/((combined['Gini Error (84%)']+combined['Gini Error (16%)'])/2),label='Gini',bins=bins,histtype='step',linewidth=2)
        ax[1,1].hist(np.abs(combined['M20']/((combined['M20 Error (84%)']+combined['M20 Error (16%)'])/2)),label='M20',bins=bins,histtype='step',linewidth=2)
        ax[0,0].set_xlabel('Concentration S/N')
        ax[0,1].set_xlabel('Asymmetry S/N')
        ax[1,0].set_xlabel('Gini S/N')
        ax[1,1].set_xlabel('M20 S/N')

        ax[0,0].set_xscale('log')
        ax[0,1].set_xscale('log')
        ax[1,0].set_xscale('log')
        ax[1,1].set_xscale('log')
    elif(type=='combined'):
        fig, ax = plt.subplots(1, figsize=(10, 8))
        ax.hist(combined['Concentration (C)']/((combined['Concentration (C) Error (84%)']+combined['Concentration (C) Error (16%)'])/2),label='Concentration',bins=bins,histtype='step',linewidth=2)
        # # # plt.xlim([-10,10])
        ax.hist(combined['Asymmetry (A)']/((combined['Asymmetry (A) Error (84%)']+combined['Asymmetry (A) Error (16%)'])/2),label='Asymmetry',bins=bins,histtype='step',linewidth=2)
        # ax[1,0].hist(asym-np.median(asym),5,alpha=0.5,label=combined)
        ax.hist(combined['Gini']/((combined['Gini Error (84%)']+combined['Gini Error (16%)'])/2),label='Gini',bins=bins,histtype='step',linewidth=2)
        # ax[1,1].hist(gini-np.median(gini),5,alpha=0.5,label=combined)
        ax.hist(np.abs(combined['M20']/((combined['M20 Error (84%)']+combined['M20 Error (16%)'])/2)),label='M20',bins=bins,histtype='step',linewidth=2)
        # m20 values are negative so need abs, others aren't
        # ax[0,1].hist(m20-np.median(m20),5,alpha=0.5,label=i[1])
        ax.legend()
        ax.set_xlabel('S/N of Measurement')
        ax.set_xscale('log')
    else:
        ## type should = compare (Default state)
        ## so do default state
        print('[---] Type for SN_ON_ERROR() function not know, defaulting to "combined" state...')
        sn_on_error()

    # ax[0,0].set_xlabel('Concentration S/N')
    # ax.legend()
    # ax[0,1].set_xlabel('Asymmetry S/N')
    # ax[1,0].set_xlabel('Gini S/N')
    # ax[1,1].set_xlabel('M20 S/N')
    # ax.set_xlabel('S/N of Measurement')
    fig.suptitle('S/N of Error on Measurements')
    # ax.set_xscale('log')
    plt.show()

# sn_on_error()
# exit()
# 
#  
def sn_on_error_compare(zoomed=False):
    ## scatter plot to compare relative error for measurements (in 2x1 plot)
    fig, ax = plt.subplots(1,2,figsize=(12,8))
    conc_x = combined['Concentration (C)']/((combined['Concentration (C) Error (84%)']+combined['Concentration (C) Error (16%)'])/2)
    asym_y = combined['Asymmetry (A)']/((combined['Asymmetry (A) Error (84%)']+combined['Asymmetry (A) Error (16%)'])/2)
    ax[0].plot(conc_x,asym_y,'b.')
    ax[0].set_xlabel('Concentration S/N')
    ax[0].set_ylabel('Asymmetry S/N')

    gini_x = combined['Gini']/((combined['Gini Error (84%)']+combined['Gini Error (16%)'])/2)
    m20_y = np.abs(combined['M20']/((combined['M20 Error (84%)']+combined['M20 Error (16%)'])/2))
    ax[1].plot(gini_x,m20_y,'b.')
    ax[1].set_xlabel('Gini S/N')
    ax[1].set_ylabel('M20 S/N')

    # ax[0].set_xscale('log')
    # ax[1].set_xscale('log')
    # ax[0].set_yscale('log')
    # ax[1].set_yscale('log')
    title = 'S/N of Error on Measurements'
    if(zoomed):
        max = 250
        ax[0].set_xlim(0,max)
        # ax[1].set_xlim(0,2000)
        # ax[1].set_xlim(0,max)
        ax[0].set_ylim(0,max)
        ax[1].set_ylim(0,max)
        ax[0].set_aspect('equal')
        # ax[1].set_aspect('equal')

        ## add the minimum tracers 
        # concentration, asymm
        ax[0].hlines(np.min(asym_y),0,max,linestyle='-.',color='r',alpha=0.3)
        ax[0].vlines(np.min(conc_x),0,max,linestyle='-.',color='r',alpha=0.3)
        ax[0].annotate(f'{np.min(asym_y):.1f}',xy={0,np.min(asym_y)},xytext=(-25,np.min(asym_y)),fontsize=12,va='center')
        ax[0].annotate(f'{np.min(conc_x):.1f}',xy={np.min(conc_x),0},xytext=(np.min(conc_x),2),fontsize=12,ha='center')

        # gini, m20
        ax[1].hlines(np.min(m20_y),0,2000,linestyle='-.',color='r',alpha=0.3)
        ax[1].vlines(np.min(gini_x),0,max,linestyle='-.',color='r',alpha=0.3)
        ax[1].annotate(f'{np.min(m20_y):.1f}',xy={0,np.min(m20_y)},xytext=(-25,np.min(m20_y)),fontsize=12,va='center')
        ax[1].annotate(f'{np.min(gini_x):.1f}',xy={np.min(gini_x),0},xytext=(np.min(gini_x),2),fontsize=12,ha='center')
        title = 'S/N of Error on Measurements (zoomed)'
    fig.suptitle(f'{title}')
    plt.tight_layout(pad=1)
    plt.show()

# sn_on_error_compare(zoomed=True)
# exit()




def source_z_compare():
    # scatter of source redshift comparisons (jades and Kirk redshift comparisons)
    fig = plt.figure(figsize=(6,5))
    plt.scatter(combined['AGN(%)'],np.abs(combined['z offset (within +- 0.05)']))
    sources = combined.query('`z offset (within +- 0.05)`>=0.2 or `z offset (within +- 0.05)`<=-0.2')
    plt.title('Z Comparison Between Catalogs')
    plt.ylabel('Z Offset')
    plt.xlabel('AGN(%)')
    plt.grid(True)
    # plt.axis('square')
    plt.ylim([0, 0.42])
    # plt.xlim([min(combined['AGN(%)']), max(combined['AGN(%)'])])
    plt.tight_layout(w_pad=1)

    plt.savefig(f'{working_directory}source_z_comparison.png',dpi=200)
    plt.show()
    plt.close()

# source_z_compare()
# exit()


agn_dom = combined.query('`AGN(%)`>=80')
comp_dom = combined.query('`AGN(%)`<80 & `AGN(%)`>=20')
sf_dom = combined.query('`AGN(%)`<20')

## different graphing groups needed
# agn_class = agn_dom.groupby('Classification')
# comp_class = comp_dom.groupby('Classification')
# sf_class = sf_dom.groupby('Classification')

agn_class = agn_dom.groupby('Main Class')
comp_class = comp_dom.groupby('Main Class')
sf_class = sf_dom.groupby('Main Class')

combined_class = combined.groupby('Main Class')
# print(combined_class.get_group('Disk')[['ID','Concentration (C)','M20','Gini','AGN(%)']])
# print(combined_class.get_group('Irregular')[['ID','Concentration (C)','M20','Gini','AGN(%)']])
# print(combined_class.get_group('Spheroid')[['ID','Concentration (C)','M20','Gini','AGN(%)']])
# exit()

def highlight_mergers(plot,x='AGN(%)',y='AGN(%)'):
    all_mergers = combined.query('`Merger (flag threshold)`=="yes"')
    plot.scatter(all_mergers[x],all_mergers[y],marker=11,label='Mergers',color='black')
    return plot

def highlight_clumps_and_arms(plot,x='AGN(%)',y='AGN(%)'):
    all_spiral_arms = combined.query('`Spiral Arms (flag)`=="yes"')
    plot.scatter(all_spiral_arms[x],all_spiral_arms[y],marker=1,label='Spiral Arms',color='purple')
    all_clumps = combined.query('`Has Clumps (flag)`=="yes"')
    plot.scatter(all_clumps[x],all_clumps[y],marker=3,label='Clumps',color='brown')
    return plot

def A_S_scatter(): # not using as of 7/16 ***
    ### scatter the asymmetry vs smoothness (clumpiness) measurements
    xlabel = 'Smoothness (S)'
    ylabel = 'Asymmetry (A)'
    savename = f'/Asymmetry_Concentration/A_S_conselice_full.png'
    xsearch=''
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))

    ## overplot the Conselice values in gray to show the different ranges
    for g in conselice_values:
        smoo = conselice_values[g][4]
        smoo_err = conselice_values[g][5]
        asymm = conselice_values[g][2]
        asymm_err = conselice_values[g][3]
        if(g=='Ellipticals'): # label it 
            plt.errorbar(smoo,asymm,xerr=smoo_err,yerr=asymm_err,fmt='.',alpha=0.5,color='darkgray',label='Conselice 2014')
        else: # dont label
            plt.errorbar(smoo,asymm,xerr=smoo_err,yerr=asymm_err,fmt='.',alpha=0.5,color='darkgray')
        plt.annotate(g,(smoo,asymm),alpha=.8)
    
    lab = ['AGN','Composite','SF']
    color_agn = ['#1f77b4','#ff7f0e','#2ca02c']
    
    # plt.scatter(agn_dom[xsearch],agn_dom[ysearch],label=lab[0],color=color_agn[0])
    # plt.scatter(comp_dom[xsearch],comp_dom[ysearch],label=lab[1],color=color_agn[1])
    # plt.scatter(sf_dom[xsearch],sf_dom[ysearch],label=lab[2],color=color_agn[2])

    b=0
    for g in [agn_class,comp_class,sf_class]:
        for c in main_colors:
            try:
                if(c=='Even Agreement'): # highlight it with * next to it or something
                    # print(f'even agreement called for {b}')
                    # id = g.get_group(c)['ID']
                    plt.annotate((f'***'),(g.get_group(c)[xsearch],g.get_group(c)[ysearch]))
                    
                    # plt.scatter(g.get_group(c)[xsearch],g.get_group(c)[ysearch],marker=main_colors[c][1],color=color_agn[b],label=(f'{lab[b]}-{c}'))
                # else:
                
                plt.scatter(g.get_group(c)[xsearch],g.get_group(c)[ysearch],marker=main_colors[c][1],color=color_agn[b],label=(f'{lab[b]}-{c}'))
                    # plt.legend()
            except:
                print(f'errrrrrrrr ----- on classification "{c}" for agn class: {b}')
                continue
            # plt.scatter(g.get_group('Mostly Spheroid')['Smoothness (S)'],g.get_group('Mostly Spheroid')['Asymmetry (A)'],marker=bin_colors['Mostly Spheroid'][1],color=color_agn[b],label=lab[b])
        b+=1
    # plt.legend()
    # plt.show()
    highlight_mergers(plt,x=xsearch,y=ysearch)
    highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)
    # plt.scatter(all_mergers['Smoothness (S)'],all_mergers['Asymmetry (A)'],marker=11,label='Mergers',color='black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.xlim([-.01,.25])
    plt.legend(loc='lower right',prop={"size":8})
    plt.tight_layout(pad=2)

    plt.title(f'{ylabel} vs {xlabel}')
    plt.gca().invert_yaxis()
    # plt.savefig(f'{working_directory}/{savename}',dpi=300)
    plt.show()
    plt.close()
    # exit()

# A_S_scatter()
# exit()

def A_C_scatter(): ###* not using as of 8/20 **
    ### scatter the asymmetry vs concentration (clumpiness) measurements
    xlabel = 'Concentration (C)'
    ylabel = 'Asymmetry (A)'
    savename = f'/Asymmetry_Concentration/A_C_conselice_full.png'
    xsearch=''
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))

    ## overplot the Conselice values in gray to show the different ranges
    # for g in conselice_values:
    #     conc = conselice_values[g][0]
    #     conc_err = conselice_values[g][1]
    #     asymm = conselice_values[g][2]
    #     asymm_err = conselice_values[g][3]
    #     if(g=='Ellipticals'): # label it 
    #         plt.errorbar(conc,asymm,xerr=conc_err,yerr=asymm_err,fmt='.',alpha=0.5,color='darkgray',label='Conselice 2014')
    #     else: # dont label
    #         plt.errorbar(conc,asymm,xerr=conc_err,yerr=asymm_err,fmt='.',alpha=0.5,color='darkgray')
    #     plt.annotate(g,(conc,asymm),alpha=.8)
    
    lab = ['AGN','Composite','SF']
    color_agn = ['#1f77b4','#ff7f0e','#2ca02c']
    
    # plt.scatter(agn_dom[xsearch],agn_dom[ysearch],label=lab[0],color=color_agn[0])
    # plt.scatter(comp_dom[xsearch],comp_dom[ysearch],label=lab[1],color=color_agn[1])
    # plt.scatter(sf_dom[xsearch],sf_dom[ysearch],label=lab[2],color=color_agn[2])

     ### plot the avg & std of each section
    b=0
    for i in [agn_dom,comp_dom,sf_dom]:
        plt.errorbar(x=(np.max(i[xsearch])+np.min(i[xsearch]))/2,y=np.median(i[ysearch]),xerr=np.std(i[xsearch]),yerr=np.std(i[ysearch]),c=color_agn[b],fmt='.',alpha=0.5)
        b+=1
    b=0
    for g in [agn_class,comp_class,sf_class]:
        for c in main_colors:
            try:
                if(c=='Even Agreement'): # highlight it with * next to it or something
                    # print(f'even agreement called for {b}')
                    # id = g.get_group(c)['ID']
                    plt.annotate((f'***'),(g.get_group(c)[xsearch],g.get_group(c)[ysearch]))
                    
                    # plt.scatter(g.get_group(c)[xsearch],g.get_group(c)[ysearch],marker=main_colors[c][1],color=color_agn[b],label=(f'{lab[b]}-{c}'))
                # else:
                
                plt.scatter(g.get_group(c)[xsearch],g.get_group(c)[ysearch],marker=main_colors[c][1],color=color_agn[b],label=(f'{lab[b]}-{c}'))
                    # plt.legend()
            except:
                print(f'errrrrrrrr ----- on classification "{c}" for agn class: {b}')
                continue
            # plt.scatter(g.get_group('Mostly Spheroid')['Smoothness (S)'],g.get_group('Mostly Spheroid')['Asymmetry (A)'],marker=bin_colors['Mostly Spheroid'][1],color=color_agn[b],label=lab[b])
        b+=1
    # plt.legend()
    # plt.show()
    highlight_mergers(plt,x=xsearch,y=ysearch)
    highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)
    # plt.scatter(all_mergers['Smoothness (S)'],all_mergers['Asymmetry (A)'],marker=11,label='Mergers',color='black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.xlim([-.01,.25])
    plt.legend(loc='lower right',prop={"size":8})
    plt.tight_layout(pad=2)

    plt.title(f'{ylabel} vs {xlabel}')
    plt.gca().invert_yaxis()
    # plt.savefig(f'{working_directory}/{savename}',dpi=300)
    plt.show()
    plt.close()
    # exit()

# A_C_scatter()
# exit()


'''
for Asymmetry standard error = 0.030654527834188688
for Concentration standard error = 0.09318416087613562
for Gini standard error = 0.01819791241163235
for M20 standard error = 0.06365918339442633
'''

# old standard error vals
standard_err_vals = {'Asymmetry (A)':0.030654527834188688,
                     'Concentration (C)':0.09318416087613562,
                     'Gini':0.01819791241163235,
                     'M20':0.06365918339442633}

## taken from messing-with-psf.py
asym_err_vals = np.array([0.0,0.0,0.0007959019999999595,0.0,-0.04180115000000001,0.0,0.0,-0.010094659000000006,0.0,-0.036242992,-0.024676754999999995,0.045978778000000026,-0.134162969,0.03602204399999999,0.0,0.0,0.0,0.0,-0.05702032000000001,-0.022280037000000003,0.0,-0.0010266109999999967,0.0,0.0,0.0,-0.04002197299999999,0.0,-0.03643905399999997,0.0,0.0031646189999999796,0.0,0.0,0.0,0.0,0.0,0.0,-0.063339852,0.0,-0.033187839,-0.06682363200000002,-0.07554662300000001,-0.012827393999999992,0.015191227000000002])
conc_err_vals = np.array([0.0,0.0,0.027574169000000204,0.0,0.11051997700000005,0.0,0.0,0.049465163000000256,0.0,0.004018027999999951,0.2278917090000001,0.039814759999999616,-0.11680277900000036,-0.11378974699999977,0.0,0.0,0.0,0.0,0.0019949949999999994,0.04489362099999994,0.0,-0.16274078000000003,0.0,0.0,0.0,-0.13341934600000016,0.0,0.09568246800000013,0.0,0.062111777000000146,0.0,0.0,0.0,0.0,0.0,0.0,0.20329700699999975,0.0,-0.06529853799999996,-0.22323508800000003,-0.15791481099999993,0.1787019970000001,0.2623501159999999])
gini_err_vals = np.array([0.0,0.0,0.01952327499999995,0.0,-0.0018238389999999938,0.0,0.0,-0.0023421839999999694,0.0,-0.01723438399999999,0.007474804999999973,-0.03887671599999998,-0.02512950900000005,-0.04304711299999997,0.0,0.0,0.0,0.0,0.023876820999999993,-0.0033005289999999965,0.0,0.0023447930000000117,0.0,0.0,0.0,-0.006689643000000023,0.0,0.010933203000000002,0.0,-0.003909011999999934,0.0,0.0,0.0,0.0,0.0,0.0,0.04591201900000008,0.0,0.0024785849999999776,-0.006363677000000012,-0.06075000600000002,0.02822096399999996,0.04588838799999995])
m20_err_vals = np.array([0.0,0.0,-0.08584868999999995,0.0,-0.23620932499999991,0.0,0.0,-0.13192586500000014,0.0,-0.016333382000000007,-0.09451803499999967,-0.11676924200000016,-0.013999436999999837,0.10628174099999987,0.0,0.0,0.0,0.0,-0.037035994000000017,-0.016739963000000024,0.0,0.09923172899999999,0.0,0.0,0.0,0.0755867939999999,0.0,-0.138279896,0.0,-0.008121830000000108,0.0,0.0,0.0,0.0,0.0,0.0,-0.1253663410000001,0.0,0.00395804100000019,-0.0007703980000000055,0.089619165,-0.10008580499999997,-0.06815098899999983])

err_vals_aligned = {'Asymmetry (A)':asym_err_vals,
                     'Concentration (C)':conc_err_vals,
                     'Gini':gini_err_vals,
                     'M20':m20_err_vals}

def make_agn_plot():
    vis_sym = {
    'Disk': '*',
    'Spheroid': 'o',
    'Irregular': '^',
    'Even Agreement': 's'
    }
    ### plot the c,a, and gini for each source by AGN fraction (on x axis)
    cols = ['green','orange','blue']
    names = ['SF','Composite','AGN']
    fig, axs = plt.subplots(3,1)
    fig.set_size_inches(10,10)
    print(f'total sources inputted to agn plot: {len(combined)}')

    def add_axis(ax_name,axis=0,lab=''):
        plot_vals = np.zeros(len(combined))
        err_vals = np.zeros((len(combined),2))
        x_vals = np.zeros(len(combined))
        plot_ind = 0
        for ind,row in combined.iterrows():
            if(row['AGN^a']>=80): c=2
            if((row['AGN^a']<80)&(row['AGN^a']>=20)): c=1
            if(row['AGN^a']<20): c=0


            plot_vals[plot_ind] = row[ax_name]
            x_vals[plot_ind] = row['AGN^a']

            ## get visual symbol here
            # print(row['Main Class'])
            # s=0
            # if(row['Main Class']=='Disk'):s=0
            # elif(row['Main Class']=='Spheroid'):s=1
            # elif(row['Main Class']=='Irregular'):s=2
            # elif(row['Main Class']=='Even Agreement'):s=3


            err_full_arr = np.zeros(100)
            err_str_arr = row[f'{ax_name} Error Full'].split('[')[1].split(']')[0].split(',')
            err_full_arr = np.array(err_str_arr)
            err_full = err_full_arr.astype(float)


            y_med = np.median(row[f'{ax_name}'])
            y_err_16 = np.abs(y_med-np.percentile(err_full,16))
            # print(y_err_16)
            y_err_84 = np.abs(y_med-np.percentile(err_full,84))
            # print(y_err_84)

            ### old version was to do standard_err_vals[ax_name]**2+(y_err_84**2), etc
            err_high = np.sqrt((standard_err_vals[ax_name]**2)+(y_err_84**2))
            err_low = np.sqrt((standard_err_vals[ax_name]**2)+(y_err_16**2))

            # err_high = np.sqrt(((err_vals_aligned[ax_name][ind])**2)+(y_err_84**2))
            # err_low = np.sqrt(((err_vals_aligned[ax_name][ind])**2)+(y_err_16**2))
            y_err = np.array(err_low,err_high)
            err_vals[plot_ind] = y_err
            axs[axis].scatter(row['AGN^a'],row[f'{ax_name}'],marker=vis_sym[row['Main Class']],color=cols[c],s=50,label=names[c],zorder=1)

            axs[axis].errorbar(row['AGN^a'],row[f'{ax_name}'],yerr=y_err,marker='.',color=cols[c],alpha=0.5,zorder=0)
            if(row['Merger (flag threshold)']=="yes"):
                axs[axis].scatter(row['AGN^a'],row[f'{ax_name}'],marker='v',label='Mergers',color='black',zorder=2)
                #  plot.scatter(all_mergers[x],all_mergers[y],marker=11,label='Mergers',color='black')
            axs[axis].set_ylabel(lab)
            axs[axis].set_xticks([])
            # locs = axs[axis].get_yticks()  # Get the current locations and labels.
            # axs[axis].set_yticks(locs,major=True,minor=True)  # Set label locations.
            plot_ind+=1
        trend_vals = np.zeros(3)
        trend_err = np.zeros(3)
        trend_x_vals = np.zeros(3)
        ic=0
        for i in [agn_dom,comp_dom,sf_dom]:
            graph_mask = np.isfinite(i[ax_name])
            # axs[axis].errorbar(x=(np.max(i['AGN(%)'][graph_mask])+np.min(i['AGN(%)'][graph_mask]))/2,y=np.median(i[ax_name][graph_mask]),yerr=np.std(i[ax_name][graph_mask]),c='gray',fmt='.',alpha=0.5)
            trend_vals[ic] = np.median(i[ax_name][graph_mask])
            trend_err[ic] = np.std(i[ax_name][graph_mask])
            trend_x_vals[ic] = (np.max(i['AGN(%)'][graph_mask])+np.min(i['AGN(%)'][graph_mask]))/2
            ic+=1
        ## plot the trend line for each
        # axs[axis].plot(trend_x_vals,trend_vals,c='red')
        # axs[axis].fill_between(trend_x_vals,trend_vals-trend_err,trend_vals+trend_err,color='red',alpha=0.5)

        # Calculating the coefficients of the best fit line 
        # graph_mask = np.isfinite(i[ax_name])
        # x_vals = combined['AGN(%)'][graph_mask]
        # coefficients = np.polyfit(x_vals, combined[ax_name][graph_mask], 1) 
        
        # plot_vals = np.abs(plot_vals)
        graph_mask = np.isfinite(plot_vals)
        x_vals = x_vals[graph_mask]
        plot_vals = plot_vals[graph_mask]
        err_vals = err_vals[graph_mask]
        arr1inds = x_vals.argsort()
        x_vals_srt = x_vals[arr1inds[::-1]]
        plot_vals_srt = plot_vals[arr1inds[::-1]]
        err_vals_srt = err_vals[arr1inds[::-1]]
        # coefficients, cov = np.polyfit(x_vals, plot_vals, 1,cov=True)
        err_vals_new = np.zeros(len(err_vals_srt))
        co = 0
        for i in err_vals_srt:
            err_vals_new[co] = np.mean(i)
            co+=1
        p, cov = np.polyfit(x_vals_srt, plot_vals_srt, 1, cov=True)
        # plt.plot(x_vals, plot_vals)
        slope = p[0]
        intercept = p[1]
        axs[axis].plot(x_vals_srt, np.polyval(p,x_vals_srt),color='black',linestyle='--')
        err = np.sqrt(np.diag(cov))


        slope_gaus = np.random.normal(slope,err[0],300)
        y_int_gaus = np.random.normal(intercept,err[1],300)
        


        # y_fit_low = (slope-err[0])*x_vals_srt + intercept-err[1]
        # y_fit_high = (slope+err[0])*x_vals_srt + intercept+err[1]
        print(f'Slope for {ax_name} = {slope} +/- {err[0]} + interc({intercept}) +/- {err[1]}')
        # axs[axis].fill_between(x_vals_srt,y_fit_low,y_fit_high,alpha=0.2,color='black')
        counter = 0
        for s in slope_gaus:
            # axs[axis].fill_between(x_vals_srt,s*x_vals_srt+(y_int_gaus[counter]-err[1]),s*x_vals_srt+(y_int_gaus[counter]+err[1]),alpha=0.1,color='gray',zorder=0)
            axs[axis].plot(x_vals_srt,s*x_vals_srt+(y_int_gaus[counter]),alpha=0.1,color='gray',zorder=0)
            counter+=1
        
        # slope = coefficients[0] 
        # intercept = coefficients[1] 
        # Generating the best fit line 
        # y_fit = slope * x_vals + intercept 
        # y_err = [slope+e * x_vals+intercept for e in err]
        # Plotting the data points and the best fit line 
        # axs[axis].plot(x_vals, y_fit, color='red',linestyle='--')
        # axs[axis].fill_between(x_vals,y_err[0],y_err[1],color='red',alpha=0.2)
        
        axs[axis].axvline(20,linestyle='--',c='lightgray',zorder=-1)
        axs[axis].axvline(80,linestyle='--',c='lightgray',zorder=-1)



    add_axis('Concentration (C)',axis=0,lab='C')
    add_axis('Asymmetry (A)',axis=1,lab='A')
    add_axis('Gini',axis=2,lab='Gini')

    # from matplotlib.lines import Line2D
    # point = Line2D([0], [0], label='Disk', marker='*', markersize=0,color='gray')
    # handles, labels = plt.gca().get_legend_handles_labels()
    # by_label = dict(zip(labels, handles))
    # axs[-1].legend(by_label.values(), by_label.keys())
    axs[-1].set_xticks(np.linspace(0,100,6))
    axs[-1].set_xlabel('AGN (%)')

    # axs[0].legend()
    fig.tight_layout(w_pad=5,h_pad=1)
    # plt.savefig('research/Analysis_graphs/CAS_vs_AGN/CA_Gini_AGN_w-o_error_only_lines.png',dpi=500)
    plt.show()

    
    return

make_agn_plot()
exit()

def testing():
    plt.scatter(0,10,marker='*',s=40,label='Disk',c='gray')
    plt.scatter(10,0,marker='o',s=40,label='Spheroid',c='gray')
    plt.scatter(10,10,marker='^',s=40,label='Irregular',c='gray')
    plt.scatter(20,10,marker='v',s=40,label='Merger',c='gray')
    plt.plot([0,0,10,10],linestyle='-',label='Lotz+2008',c='gray')
    plt.legend()
    plt.show()

# testing()
# exit()


def A_agn_scatter():
    ### scatter the given xlabel vs ylabel & include highlighting of different visual classifications
    ### give savename (include file extension (.png))
    ### xsearch & ysearch are used for the table searching (if it's different than what you're labeling it as)
    ### also highlight the agn regions (agn, composite, sf)
    xlabel = 'AGN (%)'
    ylabel = 'Asymmetry (A)'
    savename = f'/CAS_vs_AGN/A_AGN_error_bar.png'
    xsearch='AGN(%)'
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))


    # plot the conselice asymmetry values
    # plt.plot(x, y, 'k-')
    # plt.fill_between(x, y-error, y+error)

    ### plot the avg & std of each section
    for i in [agn_dom,comp_dom,sf_dom]:
        plt.errorbar(x=(np.max(i[xsearch])+np.min(i[xsearch]))/2,y=np.median(i[ysearch]),yerr=np.std(i[ysearch]),c='gray',fmt='.',alpha=0.5)
        # print(f'{np.median(i[ysearch])}+/-{np.std(i[ysearch])}')

        # ### plotting errorbar for disk and spheroid
        # disks = i.groupby('Classification').get_group('Mostly Disk')
        # spheroids = i.groupby('Classification').get_group('Mostly Spheroid')
        # plt.errorbar(x=(np.max(i[xsearch])+np.min(i[xsearch]))/3,y=np.median(disks[ysearch]),yerr=np.std(disks[ysearch]),c='blue',fmt='.',alpha=0.8)
        # plt.errorbar(x=(np.max(i[xsearch])+np.min(i[xsearch]))/4,y=np.median(spheroids[ysearch]),yerr=np.std(spheroids[ysearch]),c='orange',fmt='.',alpha=0.8)

        # print(f'Disk stats (agn,comp,sf): {np.median(disks[ysearch])}+/-{np.std(disks[ysearch])}')
        # print(f'Spheroid stats (agn,comp,sf): {np.median(spheroids[ysearch])}+/-{np.std(spheroids[ysearch])}')


    # plot the scattered values
    ind = 0
    # for g in [agn_class,comp_class,sf_class]:
    for c in main_colors:
        try:
            plt.scatter(combined_class.get_group(c)['AGN(%)'],combined_class.get_group(c)['Asymmetry (A)'],40,marker=main_colors[c][1],facecolor="none", 
                        edgecolor=main_colors[c][0],lw=1.5,label=(f'{c}'))
            
            # line for adding error bars (but they're basically invisible)
            # plt.errorbar(combined_class.get_group(c)['AGN(%)'],combined_class.get_group(c)['Asymmetry (A)'],yerr=(combined_class.get_group(c)['Asymmetry (A) Error (16%)'],combined_class.get_group(c)['Asymmetry (A) Error (84%)']),fmt='.')
            
        except:
            print(f'err with a_agn_scatter() for main_color: {c}')
            continue
    ind+=1

    # for index, row in combined.iterrows():
    #     plt.annotate(row['ID'],(row[xsearch],row[ysearch]),textcoords="offset points", xytext=(0,10), ha='center')

    # highlight the mergers (if any)
    highlight_mergers(plt,x='AGN(%)',y='Asymmetry (A)')
    # highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)

    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    # Highlight the agn regions
    y_max = plt.gca().get_ylim()[1]
    x_max = plt.gca().get_xlim()[1]

    plt.axvspan(80, 100, color='#1f77b4', alpha=0.3)
    plt.text(x=.90*x_max, y=.9*y_max, s='AGN', fontsize=12, ha='center', va='center')
    plt.axvspan(20, 80, color='#ff7f0e', alpha=0.3)
    plt.text(x=.50*x_max, y=.9*y_max, s='Composite', fontsize=12, ha='center', va='center')
    plt.axvspan(0, 20, color='#2ca02c', alpha=0.3)
    plt.text(x=.1*x_max, y=0.9*y_max, s='SF', fontsize=12, ha='center', va='center')

    plt.legend(loc='center right')
    plt.ylim((0,0.9))
    plt.title(f'{ylabel} to {xlabel} Comparison')
    plt.tight_layout(pad=2)

    # plt.savefig(f'{working_directory}/{savename}')
    # plt.gca().invert_yaxis()
    plt.show()
    plt.close()

# A_agn_scatter()
# exit()

def A_agn_histogram():
    ### histogram given xlabel vs ylabel & include highlighting of different visual classifications
    ### give savename (include file extension (.png))
    ### xsearch & ysearch are used for the table searching (if it's different than what you're labeling it as)
    ### also highlight the agn regions (agn, composite, sf)
    # xlabel = 'AGN (%)'
    ylabel = 'Count'
    xlabel = 'Asymmetry (A)'
    savename = f'/CAS_vs_AGN/A_AGN_histogram.png'
    xsearch=''
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))


    # plot the conselice asymmetry values
    # plt.plot(x, y, 'k-')
    # plt.fill_between(x, y-error, y+error)

    ### plot the avg & std of each section
    # for i in [agn_dom,comp_dom,sf_dom]:
    #     plt.errorbar(x=(np.max(i[xsearch])+np.min(i[xsearch]))/2,y=np.median(i[ysearch]),yerr=np.std(i[ysearch]),c='gray',fmt='.',alpha=0.5)
    #     print(f'{np.median(i[ysearch])}+/-{np.std(i[ysearch])}')


    lab = ['AGN','Composite','SF']
    color_agn = ['#1f77b4','#ff7f0e','#2ca02c']
    histbins = 7
    bins=np.histogram(np.hstack((agn_dom[xsearch],comp_dom[xsearch],sf_dom[xsearch])), bins=histbins)[1] #get the bin edges

    # plot the scattered values
    ind = 2
    for g in [sf_dom,comp_dom,agn_dom]:
    # for c in main_colors:
        # try:
        # plt.scatter(combined_class.get_group(c)['AGN(%)'],combined_class.get_group(c)['Asymmetry (A)'],40,marker=main_colors[c][1],facecolor="none", 
        #                 edgecolor=main_colors[c][0],lw=1.5,label=(f'{c}'))
        filling=False
        # if(ind==0):
        #     filling = True
        
        plt.hist(g[xsearch],bins=bins,label=lab[ind],color=color_agn[ind],fill=filling,histtype='step',linewidth=2)
        plt.axvline(np.median(g[xsearch]),color=color_agn[ind],linestyle="--",linewidth=2)
            # plt.legend()
        # except:
        #     print(f'err with a_agn_scatter() for main_color: {c}')
        #     continue
        ind-=1

    # for index, row in combined.iterrows():
    #     plt.annotate(row['ID'],(row[xsearch],row[ysearch]),textcoords="offset points", xytext=(0,10), ha='center')

    # highlight the mergers (if any)
    # highlight_mergers(plt,x='AGN(%)',y='Asymmetry (A)')
    # highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)

    # plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    # Highlight the agn regions
    # y_max = plt.gca().get_ylim()[1]
    # x_max = plt.gca().get_xlim()[1]

    # plt.axvspan(80, 100, color='#1f77b4', alpha=0.3)
    # plt.text(x=.90*x_max, y=.9*y_max, s='AGN', fontsize=12, ha='center', va='center')
    # plt.axvspan(20, 80, color='#ff7f0e', alpha=0.3)
    # plt.text(x=.50*x_max, y=.9*y_max, s='Composite', fontsize=12, ha='center', va='center')
    # plt.axvspan(0, 20, color='#2ca02c', alpha=0.3)
    # plt.text(x=.1*x_max, y=0.9*y_max, s='SF', fontsize=12, ha='center', va='center')

    plt.legend(loc='center right')
    # plt.ylim((0,0.9))
    plt.title(f'{xlabel} Distribution')
    plt.tight_layout(pad=1.5)

    plt.savefig(f'{working_directory}/{savename}')
    # plt.gca().invert_yaxis()
    plt.show()
    plt.close()

# A_agn_histogram()
# exit()

def C_agn_scatter():
    ### scatter the given xlabel vs ylabel & include highlighting of different visual classifications
    ### give savename (include file extension (.png))
    ### xsearch & ysearch are used for the table searching (if it's different than what you're labeling it as)
    ### also highlight the agn regions (agn, composite, sf)
    xlabel = 'AGN (%)'
    ylabel = 'Concentration (C)'
    savename = f'/CAS_vs_AGN/C_AGN_error_bar.png'
    xsearch='AGN(%)'
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))


    # plot the conselice concentration values
    color_list = [
        '#1f77b4',  # Blue
        '#ff7f0e',  # Orange
        '#2ca02c',  # Green
        '#d62728',  # Red
        '#9467bd',  # Purple
        '#8c564b',  # Brown
        '#e377c2',  # Pink
        '#17becf',  # Cyan
    ]

    # for g in conselice_values:
    #     ymin = conselice_values[g][0]-conselice_values[g][1]
    #     ymax = conselice_values[g][0]+conselice_values[g][1]
    #     c = random.choice(color_list)
    #     plt.axhspan(ymin,ymax,color=c,alpha=0.05)
    #     # plt.text(x=.60*xmax, y=(ymax-ymin)/2, s=g, fontsize=12, ha='center', va='center')
    #     plt.text(x=30, y=ymax, s=g,c=c,fontsize=12, ha='center', va='center',alpha=0.8)
    #     color_list.remove(c)


    # plt.plot(conselice_values[g][0],yerr=[g][1])
    # plt.plot(x, y, 'k-')
    # plt.fill_between(x, y-error, y+error)


    ### plot the avg & std of each section
    for i in [agn_dom,comp_dom,sf_dom]:
        graph_mask = np.isfinite(i[ysearch])
        plt.errorbar(x=(np.max(i[xsearch][graph_mask])+np.min(i[xsearch][graph_mask]))/2,y=np.median(i[ysearch][graph_mask]),yerr=np.std(i[ysearch][graph_mask]),c='gray',fmt='.',alpha=0.5)

    # plot the scattered values
    ind = 0
    # for g in [agn_class,comp_class,sf_class]:

    for c in main_colors:
        try:
            # y_err = np.array([combined_class.get_group(c)[f'{ysearch} Error (16%)'],combined_class.get_group(c)[f'{ysearch} Error (84%)']])
            # plt.errorbar(x=combined_class.get_group(c)[xsearch],y=combined_class.get_group(c)[ysearch],yerr=y_err,fmt=',')
            
            graph_mask = np.isfinite(combined_class.get_group(c)[ysearch])
            plt.scatter(combined_class.get_group(c)[xsearch][graph_mask],combined_class.get_group(c)[ysearch][graph_mask],40,marker=main_colors[c][1],facecolor="none", 
                        edgecolor=main_colors[c][0],lw=1.5,label=(f'{c}'))
            
            ## for getting errorbars (but they're basically invisible anyways)
            # plt.errorbar(combined_class.get_group(c)[xsearch][graph_mask],combined_class.get_group(c)[ysearch][graph_mask],yerr=(combined_class.get_group(c)[f'{ysearch} Error (16%)'][graph_mask],combined_class.get_group(c)[f'{ysearch} Error (84%)'][graph_mask]),fmt='.')

        except:
            print(f'err with c_agn_scatter() for main_color: {c}')
            continue
    ind+=1

    # highlight the mergers (if any)
    highlight_mergers(plt,x=xsearch,y=ysearch)


    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    # Highlight the agn regions
    y_max = plt.gca().get_ylim()[1]
    x_max = plt.gca().get_xlim()[1]

    plt.axvspan(80, 100, color='#1f77b4', alpha=0.3)
    plt.text(x=.90*x_max, y=.9*y_max, s='AGN', fontsize=12, ha='center', va='center')
    plt.axvspan(20, 80, color='#ff7f0e', alpha=0.3)
    plt.text(x=.50*x_max, y=.9*y_max, s='Composite', fontsize=12, ha='center', va='center')
    plt.axvspan(0, 20, color='#2ca02c', alpha=0.3)
    plt.text(x=.1*x_max, y=0.9*y_max, s='SF', fontsize=12, ha='center', va='center')

    plt.legend(loc='lower right')
    # plt.xlim((-.025,.125))
    plt.tight_layout(pad=1.5)
    plt.title(f'{ylabel} to {xlabel} Comparison')
    # plt.savefig(f'{working_directory}/{savename}')
    # plt.gca().invert_yaxis()
    plt.show()
    plt.close()

# C_agn_scatter()
# exit()

def S_agn_scatter(): ### not using as of 7/29 ------

    ### scatter the given xlabel vs ylabel & include highlighting of different visual classifications
    ### give savename (include file extension (.png))
    ### xsearch & ysearch are used for the table searching (if it's different than what you're labeling it as)
    ### also highlight the agn regions (agn, composite, sf)
    xlabel = 'AGN (%)'
    ylabel = 'Clumpiness (S)'
    savename = f'/CAS_vs_AGN/S_AGN_conselice_full.png'
    xsearch='AGN(%)'
    ysearch='Smoothness (S)'

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))


    # plot the conselice concentration values
    color_list = [
        '#1f77b4',  # Blue
        '#ff7f0e',  # Orange
        '#2ca02c',  # Green
        '#d62728',  # Red
        '#9467bd',  # Purple
        '#8c564b',  # Brown
        '#e377c2',  # Pink
        '#17becf',  # Cyan
    ]

    # for g in conselice_values:
    #     ymin = conselice_values[g][4]-conselice_values[g][5]
    #     ymax = conselice_values[g][4]+conselice_values[g][5]
    #     c = random.choice(color_list)
    #     plt.axhspan(ymin,ymax,color=c,alpha=0.05)
    #     # plt.text(x=.60*xmax, y=(ymax-ymin)/2, s=g, fontsize=12, ha='center', va='center')
    #     plt.text(x=30, y=ymin, s=g,c=c,fontsize=12, ha='center', va='center',alpha=0.8)
    #     color_list.remove(c)


    # plt.plot(conselice_values[g][0],yerr=[g][1])
    # plt.plot(x, y, 'k-')
    # plt.fill_between(x, y-error, y+error)

    # plot the scattered values
    ind = 0
    # for g in [agn_class,comp_class,sf_class]:
    for c in bin_colors:
        try:
            plt.scatter(combined_class.get_group(c)[xsearch],combined_class.get_group(c)[ysearch],marker=bin_colors[c][1],label=(f'{c}'))
            # plt.legend()
        except:
            continue
    ind+=1

    # highlight the mergers (if any)
    highlight_mergers(plt,x=xsearch,y=ysearch)


    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    # Highlight the agn regions
    # plt.ylim((-.1,.4))
    y_max = plt.gca().get_ylim()[1]
    x_max = plt.gca().get_xlim()[1]

    plt.axvspan(80, 100, color='lightgray', alpha=0.3)
    plt.text(x=.90*x_max, y=.9*y_max, s='AGN', fontsize=12, ha='center', va='center')
    plt.axvspan(20, 80, color='lightgray', alpha=0.3)
    plt.text(x=.50*x_max, y=.9*y_max, s='Composite', fontsize=12, ha='center', va='center')
    plt.axvspan(0, 20, color='lightgray', alpha=0.3)
    plt.text(x=.1*x_max, y=0.9*y_max, s='SF', fontsize=12, ha='center', va='center')

    plt.legend(loc='best',prop={'size':8})
    plt.tight_layout(pad=2)
    plt.title(f'{nir_filter.upper()} {ylabel} to {xlabel} Comparison')
    plt.savefig(f'{working_directory}/{savename}')
    plt.gca().invert_yaxis()
    plt.show()
    plt.close()

def gini_agn_scatter():
    ### scatter the given xlabel vs ylabel & include highlighting of different visual classifications
    ### give savename (include file extension (.png))
    ### xsearch & ysearch are used for the table searching (if it's different than what you're labeling it as)
    ### also highlight the agn regions (agn, composite, sf)
    xlabel = 'AGN(%)'
    # xlabel = 'Asymmetry (A)'
    ylabel = 'Gini'
    savename = f'/CAS_vs_AGN/gini_AGN_full.png'
    xsearch=''
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))


    ## plot the avg & std of each section
    for i in [agn_dom,comp_dom,sf_dom]:
        graph_mask = np.isfinite(i[ysearch])
        plt.errorbar(x=(np.max(i[xsearch][graph_mask])+np.min(i[xsearch][graph_mask]))/2,y=np.median(i[ysearch][graph_mask]),yerr=np.std(i[ysearch][graph_mask]),c='gray',fmt='.',alpha=0.5)

    # plot the scattered values
    ind = 0
    # for g in [agn_class,comp_class,sf_class]:
    for c in main_colors:
        try:
            graph_mask = np.isfinite(combined_class.get_group(c)[ysearch])
            plt.scatter(combined_class.get_group(c)[xsearch][graph_mask],combined_class.get_group(c)[ysearch][graph_mask],40,marker=main_colors[c][1],facecolor="none", 
                        edgecolor=main_colors[c][0],lw=1.5,label=(f'{c}'))
            
            ## for getting errorbars
            plt.errorbar(combined_class.get_group(c)[xsearch][graph_mask],combined_class.get_group(c)[ysearch][graph_mask],yerr=(combined_class.get_group(c)[f'{ysearch} Error (16%)'][graph_mask],combined_class.get_group(c)[f'{ysearch} Error (84%)'][graph_mask]),fmt='.')
        except:
            print(f'err with gini_agn_scatter() for main_color: {c}')
            continue
    ind+=1

    # highlight the mergers (if any)
    highlight_mergers(plt,x=xsearch,y=ysearch)
    # highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)


    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    # Highlight the agn regions
    y_max = plt.gca().get_ylim()[1]
    x_max = plt.gca().get_xlim()[1]

    plt.axvspan(80, 100, color='#1f77b4', alpha=0.3)
    plt.text(x=.90*100, y=.9*y_max, s='AGN', fontsize=12, ha='center', va='center')
    plt.axvspan(20, 80, color='#ff7f0e', alpha=0.3)
    plt.text(x=.50*100, y=.9*y_max, s='Composite', fontsize=12, ha='center', va='center')
    plt.axvspan(0, 20, color='#2ca02c', alpha=0.3)
    plt.text(x=.1*100, y=0.9*y_max, s='SF', fontsize=12, ha='center', va='center')

    plt.legend(loc='center right')
    plt.tight_layout(pad=1.5)

    # plt.xlim((-.025,.125))
    plt.title(f'{ylabel} to {xlabel} Comparison')
    plt.savefig(f'{working_directory}/{savename}')
    # plt.gca().invert_yaxis()
    plt.show()
    plt.close()

# gini_agn_scatter()
# exit()

def gini_agn_histogram():
    ### histogram given xlabel vs ylabel & include highlighting of different visual classifications
    ### give savename (include file extension (.png))
    ### xsearch & ysearch are used for the table searching (if it's different than what you're labeling it as)
    ### also highlight the agn regions (agn, composite, sf)
    # xlabel = 'AGN (%)'
    ylabel = 'Count'
    xlabel = 'Gini'
    savename = f'/CAS_vs_AGN/gini_AGN_histogram.png'
    xsearch=''
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))


    # plot the conselice asymmetry values
    # plt.plot(x, y, 'k-')
    # plt.fill_between(x, y-error, y+error)

    ### plot the avg & std of each section
    # for i in [agn_dom,comp_dom,sf_dom]:
    #     plt.errorbar(x=(np.max(i[xsearch])+np.min(i[xsearch]))/2,y=np.median(i[ysearch]),yerr=np.std(i[ysearch]),c='gray',fmt='.',alpha=0.5)
    #     print(f'{np.median(i[ysearch])}+/-{np.std(i[ysearch])}')


    lab = ['AGN','Composite','SF']
    color_agn = ['#1f77b4','#ff7f0e','#2ca02c']
    histbins = 7
    graph_masks = [np.isfinite(g[xsearch]) for g in [agn_dom,comp_dom,sf_dom]]
    bins=np.histogram(np.hstack((agn_dom[xsearch][graph_masks[0]],comp_dom[xsearch][graph_masks[1]],sf_dom[xsearch][graph_masks[2]])), bins=histbins)[1] #get the bin edges
    # bins=[2,5,7]
    # plot the scattered values
    ind = 2
    for g in [sf_dom,comp_dom,agn_dom]:
    # for c in main_colors:
        # try:
        # plt.scatter(combined_class.get_group(c)['AGN(%)'],combined_class.get_group(c)['Asymmetry (A)'],40,marker=main_colors[c][1],facecolor="none", 
        #                 edgecolor=main_colors[c][0],lw=1.5,label=(f'{c}'))
        filling=False
        wdth=2+ind
        if(ind==0):
            wdth=2
        
        graph_mask = np.isfinite(g[xsearch])
        plt.hist(g[xsearch][graph_mask],bins=bins,label=lab[ind],color=color_agn[ind],fill=filling,histtype='step',linewidth=2)
            # plt.legend()
        plt.axvline(np.median(g[xsearch][graph_mask]),color=color_agn[ind],linestyle='--',linewidth=2)
        # except:
        #     print(f'err with a_agn_scatter() for main_color: {c}')
        #     continue
        ind-=1


    # for index, row in combined.iterrows():
    #     plt.annotate(row['ID'],(row[xsearch],row[ysearch]),textcoords="offset points", xytext=(0,10), ha='center')

    # highlight the mergers (if any)
    # highlight_mergers(plt,x='AGN(%)',y='Asymmetry (A)')
    # highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)

    # plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    # Highlight the agn regions
    # y_max = plt.gca().get_ylim()[1]
    # x_max = plt.gca().get_xlim()[1]

    # plt.axvspan(80, 100, color='#1f77b4', alpha=0.3)
    # plt.text(x=.90*x_max, y=.9*y_max, s='AGN', fontsize=12, ha='center', va='center')
    # plt.axvspan(20, 80, color='#ff7f0e', alpha=0.3)
    # plt.text(x=.50*x_max, y=.9*y_max, s='Composite', fontsize=12, ha='center', va='center')
    # plt.axvspan(0, 20, color='#2ca02c', alpha=0.3)
    # plt.text(x=.1*x_max, y=0.9*y_max, s='SF', fontsize=12, ha='center', va='center')

    plt.legend(loc='upper right')
    # plt.ylim((0,0.9))
    plt.title(f'{xlabel} Distribution')
    plt.tight_layout(pad=1.5)

    plt.savefig(f'{working_directory}/{savename}')
    # plt.gca().invert_yaxis()
    plt.show()
    plt.close()

# gini_agn_histogram()
# exit()


def blank_agn_plot(): # template for making the cutout plot of type and agn
    plt.xlabel('AGN (%)')
    # plt.ylabel('Classification',labelpad=30)
    plt.xlim([0,100])
    plt.yticks(visible=False)
    plt.yticks([0,1,2,3])
    plt.ylim([0,3])
    plt.tight_layout(pad=1.5)
    plt.title('Galaxy Types to AGN (%)')
    plt.savefig(f'{working_directory}cutout-agn-plot.png',dpi=300)
    plt.show()
# blank_agn_plot()
# exit()

def gini_m20_scatter_new():
    vis_sym = {
    'Disk': '*',
    'Spheroid': 'o',
    'Irregular': '^',
    'Even Agreement': 's'
    }
    ### plot the c,a, and gini for each source by AGN fraction (on x axis)
    cols = ['green','orange','blue']
    names = ['SF','Composite','AGN']
    fig, axs = plt.subplots()
    fig.set_size_inches(10,8)
    print(f'total sources inputted to agn plot: {len(combined)}')

    def add_axis(ax_name,xsearch='',axis=0,lab=''):
        gen_err=False
        for ind,row in combined.iterrows():
            if(row['AGN^a']>=80): c=2
            if((row['AGN^a']<80)&(row['AGN^a']>=20)): c=1
            if(row['AGN^a']<20): c=0

            ## get visual symbol here
            # print(row['Main Class'])
            # s=0
            # if(row['Main Class']=='Disk'):s=0
            # elif(row['Main Class']=='Spheroid'):s=1
            # elif(row['Main Class']=='Irregular'):s=2
            # elif(row['Main Class']=='Even Agreement'):s=3

            def get_err(val):
                err_full_arr = np.zeros(100)
                err_str_arr = row[f'{val} Error Full'].split('[')[1].split(']')[0].split(',')
                err_full_arr = np.array(err_str_arr)
                err_full = err_full_arr.astype(float)


                y_med = np.median(row[f'{val}'])
                y_err_16 = np.abs(y_med-np.percentile(err_full,16))
                # print(y_err_16)
                y_err_84 = np.abs(y_med-np.percentile(err_full,84))
                # print(y_err_84)
                err_high = np.sqrt((standard_err_vals[ax_name]**2)+(y_err_84**2))
                err_low = np.sqrt((standard_err_vals[ax_name]**2)+(y_err_16**2))
                y_err = np.array(err_low,err_high)
                return y_err
            y_err = get_err(ax_name)
            x_err = get_err(xsearch)
            
            axs.scatter(row[f'{xsearch}'],row[f'{ax_name}'],marker=vis_sym[row['Main Class']],color=cols[c],s=50,label=names[c],zorder=1)
            if(not gen_err): # make general errorbar if none already
                axs.errorbar(-1,0.7,xerr=x_err,yerr=y_err,fmt='-',color='black',zorder=0,capsize=6)
                gen_err=True
            if(row['Merger (flag threshold)']=="yes"):
                axs.scatter(row[xsearch],row[f'{ax_name}'],marker='v',label='Mergers',color='gray',zorder=2)
                #  plot.scatter(all_mergers[x],all_mergers[y],marker=11,label='Mergers',color='black')
            axs.set_ylabel(lab)
            # axs.set_xticks([])
            # locs = axs[axis].get_yticks()  # Get the current locations and labels.
            # axs[axis].set_yticks(locs,major=True,minor=True)  # Set label locations.
        agn_pop = combined[combined['AGN(%)']>=80]
        comp_pop = combined[(combined['AGN(%)']>=20)&(combined['AGN(%)']<80)]
        sf_pop = combined[combined['AGN(%)']<20]
        b=0
        # for i in [agn_dom,comp_dom,sf_dom]:
        for i in [sf_pop,comp_pop,agn_pop]:
            xs = i[xsearch].dropna(how='any')
            ys = i[ax_name].dropna(how='any')
            axs.errorbar(x=np.median(xs),y=np.median(ys),xerr=np.std(xs),yerr=np.std(ys),c=cols[b],fmt='s',alpha=0.5)
            # graph_mask = np.isfinite(i[ax_name])
            # x_mask = np.isfinite(i[xsearch])
            # print(i[ax_name][graph_mask])
            # print(i[xsearch])
            # axs.errorbar(x=np.mean(i[xsearch]),y=np.mean(i[ax_name]),xerr=np.std(i[xsearch]),yerr=np.std(i[ax_name]),c=cols[b],fmt='s',alpha=0.5)
            b+=1

        # print(np.mean(agn_pop[xsearch]))
        # print(agn_pop['Gini'])

        # axs.errorbar(x=np.median(i[xsearch][x_mask]),y=np.median(i[ax_name][graph_mask]),xerr=np.std(i[xsearch][graph_mask]),yerr=np.std(i[ax_name][graph_mask]),c=cols[b],fmt='s',alpha=0.5)



    add_axis('Gini',xsearch='M20',lab='Gini')

    ## show the Lotz et al. 2008 merger line
    x = [0,-3.0]
    y = [.33,.75]
    plt.plot(x, y,color='darkgray',label='Lotz+2008')
    plt.annotate('Merger',(-1.2,0.65),color='darkgray')

    ## show E/S0/Sa & Sb/Sc/Ir separator line (Lotz et al. 2008)
    x = [-1.7,-3.0]
    y = [0.568,0.38]
    plt.plot(x, y,color='darkgray')
    plt.annotate('E/S0/Sa',(-2.5,0.55),color='darkgray')
    plt.annotate('Sb/Sc/Ir',(-1.2,0.42),color='darkgray')

    plt.xlabel('M20')
    plt.tight_layout(pad=1)

    plt.xlim([-0.8,-2.8])
    plt.ylim([0.38,0.73])
    plt.savefig('research/Analysis_graphs/gini_m20/gini_m20_no_legend.png',dpi=500)
    plt.show()



gini_m20_scatter_new()
exit()

def gini_m20_scatter():
    ### scatter the asymmetry vs smoothness (clumpiness) measurements
    xlabel = '$M_{20}$'
    ylabel = 'Gini'
    savename = f'/gini_m20/gini_m20_full.png'
    xsearch='M20'
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))

    
    lab = ['AGN','Composite','SF']
    color_agn = ['blue','orange','green']
    
    
    ## show the Lotz et al. 2008 merger line
    x = [0,-3.0]
    y = [.33,.75]
    plt.plot(x, y,color='darkgray',label='Lotz+2008')
    plt.annotate('Merger',(-1.2,0.65),color='darkgray')

    ## show E/S0/Sa & Sb/Sc/Ir separator line (Lotz et al. 2008)
    x = [-1.7,-3.0]
    y = [0.568,0.38]
    plt.plot(x, y,color='darkgray')
    plt.annotate('E/S0/Sa',(-2.5,0.55),color='darkgray')
    plt.annotate('Sb/Sc/Ir',(-1.2,0.42),color='darkgray')
    # plt.scatter(sf_dom[xsearch],agn_dom[ysearch],label=lab[0],color=color_agn[0])
    # plt.scatter(comp_dom[xsearch],comp_dom[ysearch],label=lab[1],color=color_agn[1])
    # plt.scatter(sf_dom[xsearch],sf_dom[ysearch],label=lab[2],color=color_agn[2])

    ### plot the avg & std of each section
    b=0
    for i in [agn_dom,comp_dom,sf_dom]:
        graph_mask = (np.isfinite(i[ysearch])) & (np.isfinite(i[xsearch]))
        # x_mask = np.isfinite(i[xsearch])
        plt.errorbar(x=np.median(i[xsearch][graph_mask]),y=np.median(i[ysearch][graph_mask]),xerr=np.std(i[xsearch][graph_mask]),yerr=np.std(i[ysearch][graph_mask]),c=color_agn[b],fmt='.',alpha=0.5)
        b+=1


    b=0
    for g in [agn_class,comp_class,sf_class]:
        for c in main_colors:
            try:
                graph_mask = (np.isfinite(g.get_group(c)[ysearch])) & (np.isfinite(g.get_group(c)[xsearch]))

                plt.scatter(g.get_group(c)[xsearch][graph_mask],g.get_group(c)[ysearch][graph_mask],40,marker=main_colors[c][1], edgecolor=color_agn[b],lw=1.5,label=(f'{lab[b]}-{c}'))

                ## for getting errorbars
                # plt.errorbar(g.get_group(c)[xsearch][graph_mask],g.get_group(c)[ysearch][graph_mask],xerr=(g.get_group(c)[f'{xsearch} Error (16%)'][graph_mask],g.get_group(c)[f'{xsearch} Error (84%)'][graph_mask]),
                #                 yerr=(g.get_group(c)[f'{ysearch} Error (16%)'][graph_mask],g.get_group(c)[f'{ysearch} Error (84%)'][graph_mask]),fmt='.',color=color_agn[b])
            except:
                continue
        b+=1
    # plt.legend()
    # plt.show()

    highlight_mergers(plt,x=xsearch,y=ysearch)
    # highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim([-0.5,-2.8])
    # plt.ylim([0.3,0.7])
    plt.legend(loc='best',prop={"size":8.5})
    plt.tight_layout(pad=2)

    plt.title(f'{ylabel} vs {xlabel}')
    # plt.gca().invert_xaxis()
    # plt.savefig(f'{working_directory}/{savename}',dpi=300)
    plt.show()
    plt.close()

# gini_m20_scatter()
# exit()

def gini_C_scatter():
    ### scatter the asymmetry vs gini
    xlabel = 'Concentration (C)'
    ylabel = 'Gini'
    savename = f'/gini_concentration_full.png'
    xsearch=''
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))

    
    lab = ['AGN','Composite','SF']
    color_agn = ['#1f77b4','#ff7f0e','#2ca02c']
    
    ### plot the avg & std of each section
    b=0
    for i in [agn_dom,comp_dom,sf_dom]:
        graph_mask = (np.isfinite(i[ysearch])) & (np.isfinite(i[xsearch]))

        plt.errorbar(x=np.median(i[xsearch][graph_mask]),y=np.median(i[ysearch][graph_mask]),xerr=np.std(i[xsearch][graph_mask]),yerr=np.std(i[ysearch][graph_mask]),c=color_agn[b],fmt='.',alpha=0.5)
        b+=1


    b=0
    for g in [agn_class,comp_class,sf_class]:
        for c in main_colors:
            try:
                graph_mask = (np.isfinite(g.get_group(c)[ysearch])) & (np.isfinite(g.get_group(c)[xsearch]))

                plt.scatter(g.get_group(c)[xsearch][graph_mask],g.get_group(c)[ysearch][graph_mask],40,marker=main_colors[c][1],facecolor="none", edgecolor=color_agn[b],lw=1.5,label=(f'{lab[b]}-{c}'))
                ## for getting errorbars
                plt.errorbar(g.get_group(c)[xsearch][graph_mask],g.get_group(c)[ysearch][graph_mask],xerr=(g.get_group(c)[f'{xsearch} Error (16%)'][graph_mask],g.get_group(c)[f'{xsearch} Error (84%)'][graph_mask]),
                                yerr=(g.get_group(c)[f'{ysearch} Error (16%)'][graph_mask],g.get_group(c)[f'{ysearch} Error (84%)'][graph_mask]),fmt='.',color=color_agn[b])

            except:
                continue
        b+=1
    # plt.legend()
    # plt.show()

    highlight_mergers(plt,x=xsearch,y=ysearch)
    # highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.xlim([-0.5,-2.8])
    # plt.ylim([0.3,0.7])
    plt.legend(loc='best',prop={"size":8.5})
    plt.tight_layout(pad=2)

    plt.title(f'{ylabel} vs {xlabel}')
    # plt.gca().invert_xaxis()
    # plt.savefig(f'{working_directory}/{savename}',dpi=300)
    plt.show()
    plt.close()

# gini_C_scatter()
# exit()

def gini_A_scatter(): # haven't updated with error bars (10/26/24)
    ### scatter the asymmetry vs gini
    xlabel = 'Asymmetry (A)'
    ylabel = 'Gini'
    savename = f'/gini_asymmetry_full.png'
    xsearch=''
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))

    
    lab = ['AGN','Composite','SF']
    color_agn = ['#1f77b4','#ff7f0e','#2ca02c']
    
    ### plot the avg & std of each section
    # b=0
    # for i in [agn_dom,comp_dom,sf_dom]:
    #     plt.errorbar(x=(np.max(i[xsearch])+np.min(i[xsearch]))/2,y=np.median(i[ysearch]),xerr=np.std(i[xsearch]),yerr=np.std(i[ysearch]),c=color_agn[b],fmt='.',alpha=0.5)
    #     b+=1


    b=0
    for g in [agn_class,comp_class,sf_class]:
        for c in main_colors:
            try:
                plt.scatter(g.get_group(c)[xsearch],g.get_group(c)[ysearch],40,marker=main_colors[c][1],facecolor="none", edgecolor=color_agn[b],lw=1.5,label=(f'{lab[b]}-{c}'))
                # plt.legend()
            except:
                continue
        b+=1
    # plt.legend()
    # plt.show()

    highlight_mergers(plt,x=xsearch,y=ysearch)
    # highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.xlim([-0.5,-2.8])
    # plt.ylim([0.3,0.7])
    plt.legend(loc='best',prop={"size":8.5})
    plt.tight_layout(pad=2)

    plt.title(f'{ylabel} vs {xlabel}')
    # plt.gca().invert_xaxis()
    plt.savefig(f'{working_directory}/{savename}',dpi=300)
    plt.show()
    plt.close()

# gini_A_scatter()
# exit()

def m20_A_scatter():
    ### scatter the m20 & concentration
    xlabel = 'Asymmetry (A)'
    ylabel = '$M_{20}$'
    savename = f'/m20_asymmetry_error_bar.png'
    xsearch=''
    ysearch='M20'

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))

    
    lab = ['AGN','Composite','SF']
    color_agn = ['#1f77b4','#ff7f0e','#2ca02c']
    
    ### plot the avg & std of each section
    b=0
    for i in [agn_dom,comp_dom,sf_dom]:
        graph_mask = (np.isfinite(i[ysearch])) & (np.isfinite(i[xsearch]))

        plt.errorbar(x=np.median(i[xsearch][graph_mask]),y=np.median(i[ysearch][graph_mask]),xerr=np.std(i[xsearch][graph_mask]),yerr=np.std(i[ysearch][graph_mask]),c=color_agn[b],fmt='.',alpha=0.5)
        b+=1


    b=0
    for g in [agn_class,comp_class,sf_class]:
        for c in main_colors:
            try:
                graph_mask = (np.isfinite(g.get_group(c)[ysearch])) & (np.isfinite(g.get_group(c)[xsearch]))

                plt.scatter(g.get_group(c)[xsearch][graph_mask],g.get_group(c)[ysearch][graph_mask],40,marker=main_colors[c][1],facecolor="none", edgecolor=color_agn[b],lw=1.5,label=(f'{lab[b]}-{c}'))
                ## for getting errorbars
                plt.errorbar(g.get_group(c)[xsearch][graph_mask],g.get_group(c)[ysearch][graph_mask],xerr=(g.get_group(c)[f'{xsearch} Error (16%)'][graph_mask],g.get_group(c)[f'{xsearch} Error (84%)'][graph_mask]),
                                yerr=(g.get_group(c)[f'{ysearch} Error (16%)'][graph_mask],g.get_group(c)[f'{ysearch} Error (84%)'][graph_mask]),fmt='.',color=color_agn[b])

            except:
                continue
        b+=1
    # plt.legend()
    # plt.show()

    highlight_mergers(plt,x=xsearch,y=ysearch)
    # highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.xlim([-0.5,-2.8])
    # plt.ylim([0.3,0.7])
    plt.legend(loc='best',prop={"size":8.5})
    plt.tight_layout(pad=2)

    plt.title(f'{ylabel} vs {xlabel}')
    # plt.gca().invert_xaxis()
    # plt.savefig(f'{working_directory}/{savename}',dpi=300)
    plt.show()
    plt.close()

# m20_A_scatter()
# exit()


#### older graphs that haven't been update (7/29)
# type_z_hist()
# testing_agn_frac_z()
# agn_frac_merger_scatter()
# agn_frac_z()
# A_S_scatter() # not using this one as of 7/29
# S_agn_scatter() # not using this one as of 7/29
# A_C_scatter() # not using this one as of 7/29

######################################################
#### used graphs as of (7/29)
# agn_frac_hist()
# agn_frac_merger_hist()
# compare_subsample()

## non parametric measurements
A_agn_histogram()
A_agn_scatter()
C_agn_scatter()
gini_agn_histogram()
gini_agn_scatter()
gini_C_scatter()
m20_A_scatter()
gini_m20_scatter()





def plot_AGN_statmorph(xlabel,ylabel,savename,xsearch='',ysearch=''):
    ### template function for creating the statmorph comparison graphs
    # scatter the given xlabel vs ylabel & include highlighting of different visual classifications
    ### give savename (include file extension (.png))
    ### xsearch & ysearch are used for the table searching (if it's different than what you're labeling it as)
    ### also highlight the agn regions (agn, composite, sf)
    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch

    ind = 0
    for i in classes:
        if(not i.empty):
            plt.scatter(i[xsearch],stats_table[ind][ysearch],marker=symbols[ind],label=labels[ind])
        ind+=1

    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    # Highlight the agn regions
    y_max = plt.gca().get_ylim()[1]
    x_max = plt.gca().get_xlim()[1]

    plt.axvspan(80, 100, color='lightgray', alpha=0.3)
    plt.text(x=.90*x_max, y=.9*y_max, s='AGN', fontsize=12, ha='center', va='center')
    plt.axvspan(20, 80, color='lightgray', alpha=0.3)
    plt.text(x=.50*x_max, y=.9*y_max, s='Composite', fontsize=12, ha='center', va='center')
    plt.axvspan(0, 20, color='lightgray', alpha=0.3)
    plt.text(x=.1*x_max, y=0.9*y_max, s='SF', fontsize=12, ha='center', va='center')

    plt.legend(loc='center right')
    plt.title(f'{nir_filter.upper()} {ylabel} to {xlabel} Comparison')
    plt.savefig(f'{working_directory}/{savename}')
    plt.show()
    plt.close()

# plot_AGN_statmorph('AGN(%)','Gini',f'{nir_filter}AGN_gini_test.png')





disclaimer(file)
disclaimer(stat_measures)