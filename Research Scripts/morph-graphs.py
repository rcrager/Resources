import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


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
small_size = 11
large_size = 15
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
stat_measures = f"research/statmorph_output/grizli/grizli-error-bar-measurements-new-testing-9-14.tsv"
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

    plt.scatter(table['z'],table['Obs WV (um)']/(1+table['z']),marker='o')
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
    plt.axvspan(3.0,4.0,facecolor='brown',alpha=0.3) # f444w
    plt.axvspan(2.4,3.0,facecolor='purple',alpha=0.3) # f356w
    plt.axvspan(1.5,2.4,facecolor='red',alpha=0.3) # f277w (this ones kind of large and off)
    plt.axvspan(0.75,1.5,facecolor='green',alpha=0.3) # f200w
    plt.axvspan(0,0.75,facecolor='orange',alpha=0.3) # f150w
    
    plt.annotate('f444w',(3.5,4.4),horizontalalignment='center',verticalalignment='top')
    plt.annotate('f356w',(2.65,4.4),horizontalalignment='center',verticalalignment='top')
    plt.annotate('f277w',(2.0,4.4),horizontalalignment='center',verticalalignment='top')
    plt.annotate('f200w',(1.0,4.4),horizontalalignment='center',verticalalignment='top')
    plt.annotate('f150w',(.4,4.4),horizontalalignment='center',verticalalignment='top')

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



def agn_frac_merger_hist():
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
        ax[0,0].plot(i[1]['S/N'],(i[1]['Concentration (C) Error (16%)']+i[1]['Concentration (C) Error (84%)']),'.',label=f'{i[0]}')

        # asymmetry
        ax[0,1].plot(i[1]['S/N'],(i[1]['Asymmetry (A) Error (16%)']+i[1]['Asymmetry (A) Error (84%)']),'.',label=f'{i[0]}')

        # gini
        ax[1,0].plot(i[1]['S/N'],(i[1]['Gini Error (16%)']+i[1]['Gini Error (84%)']),'.',label=f'{i[0]}')

        # m20
        ax[1,1].plot(i[1]['S/N'],(i[1]['M20 Error (16%)']+i[1]['M20 Error (84%)']),'.',label=f'{i[0]}')

    # ax[0,0] = concentration
    ax[0,0].legend()
    # ax[0,0].set_title('Concentration error (*100) vs S/N by Filter')
    ax[0,0].set_ylabel('Concentration (C) Error')
    ax[0,0].set_xlabel('S/N (avg per pixel)')
    # plt.xlim([0,100])

    # ax[0,1] = asymmetry
    # ax[0,1].legend()
    # ax[0,0].title('Asymmetry error (*100) vs S/N by Filter')
    ax[0,1].set_ylabel('Asymmetry (A) Error')
    # ax[0,1].set_xlabel('S/N (avg per pixel)')

    # ax[1,0] = gini
    # ax[1,0].legend()
    ax[1,0].set_ylabel('Gini Error')
    # ax[1,0].set_xlabel('S/N (avg per pixel)')

    # ax[0,1] = asymmetry
    # ax[1,1].legend()
    ax[1,1].set_ylabel('M20 Error')
    # ax[1,1].set_xlabel('S/N (avg per pixel)')

    # ax[0,0].set_xlim([0,100])
    # ax[1,0].set_xlim([0,100])
    # ax[1,1].set_xlim([0,100])
    # ax[0,1].set_xlim([0,100])
    # plt.tight_layout()
    fig.tight_layout(pad=2)
    fig.suptitle('Relative Error (16th + 84th) vs S/N (per Filter)',fontsize=16)
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



    ## for another source 
    # Concentration for first source = [2.5130088831624384, 1.9120295449878268, 1.272858166928092, 1.812445149498708, 2.5964915422459915, 2.7750766813192818, 4.166558657443729, 1.959018683368927, 0.70474726018198, 1.650737922469481, 2.0612876482745355, 1.9817308485293006, 1.5856142498367736, 1.4382910923349916, 1.2004573109036512, 0.8784050656390991, 1.6880115613169022, 2.7033561551763237, 1.9556036236418808, 1.4560020678371846, 1.6396375526777134, 2.591053018109198, 1.5128345776646854, 1.4831577445387483, 1.1246154784488158, 1.524767268609217, 2.303498519122887, 1.3531585543701898, 1.7805909103895519, 2.1971840235481706, 0.0, 1.4235415674988583, 2.2221480542306438, 2.7471533089342817, 2.0820167290008, 7.031716483845184, 1.311699178921554, 2.037879701532675, 2.3102169061509983, 2.1549040881876707, 1.4360284764081173, 2.819899010013706, 2.297243976227889, 2.4019602640057216, 1.4186986000378812, 0.7847751773432465, 1.5628669307146834, 1.9963051195136219, 1.7814357073450668, 0.7828656214272189, 1.0630526322369274, 0.35198790042793743, 2.7853730241872263, 2.980010692742376, 1.9703780874525374, 2.2168760392904705, 1.8556563888551902, 1.4943726368386576, 1.5051500428558606, 1.4163922997431309, 2.3531820477847902, 1.7044979401806295, 1.2756300210460514, 1.049577332262969, 1.795168122533366, 1.6756061121579313, 2.919785325922078, 2.901732120663346, 2.1196964913352963, 1.5234071595283247, 2.275130606735121, 1.0672307080703827, 2.3043148869101824, 2.1853525914972867, 1.6637376406806228, 2.4455716436007475, 1.0492902435860914, 2.3044491169913033, 1.1799756120272142, 1.5223705748292817, 1.9875032905611545, 0.6701077063489541, 2.1151227603120777, 1.5051502103698036, 2.904957911475723, 2.7717607175049253, 2.5480769181522605, 0.7814443502590368, 1.7208769009348666, 3.1388295189419764, 2.4450972720431423, 0.0, 3.590356072227852, 1.9284685018362893, 1.2245349714780387, 2.659871499456801, 1.108280804553358, 0.8964047613845059, 2.882238640105542, 2.807077063742632]
    # Asymmetry for first source = [-0.62020014714609, -0.6729995636668931, -0.421289477251502, -0.32715671436036214, -0.39623954543192486, -0.3721338802017393, -0.47180441811190776, -0.41913465724580534, -0.12585870163436222, -0.3361524165183649, -0.36946559592636263, -0.48999147432492823, -0.38734883508325973, -0.39151266317635625, -0.3019913846429439, -0.5032364744759075, -0.4130609186239919, -0.41801332121834656, -0.3687952142827253, -0.3788144947857085, -0.19444579243422844, -0.6607964278283612, -0.5359511815581108, -0.4602619052791239, -0.4012788591305001, -0.30532568844638414, -0.903182290080471, -0.36470160415461, -0.4861888684861692, -0.41596022594336085, 100.0, -0.3881870649920057, -0.6618018408726954, -0.49011722449487755, -0.3835666709316191, -0.6421792905343764, -0.4963625143152089, -0.6573179521947415, -0.559330089451383, -0.7313264231731978, -0.8663599212584125, -0.7970324468909462, -0.4835970516539534, -0.4612969745181021, 0.2951903643460362, -0.39062556873659127, -0.7319915792464159, -0.6394895412111706, -0.3973399427037772, -0.41335322379499373, -0.6207150823546141, -1.3421131329546898, -0.4862189179237855, -0.6063190754009423, -0.44438225188299446, -0.5077730555811103, -0.309389881904017, -1.073257724158147, 0.008831187507142282, -0.42441144581254436, -0.36044730234591965, -0.3236579109240068, -0.3676349188076344, -0.4406694907820676, -0.4612523542340657, -0.7511294252188623, -0.5684834490588677, -0.6514297287496229, -0.37021380019745426, -0.35368649020501985, -0.6729963497788451, -0.4484402537696136, -0.35128418906934705, -0.6777723954434722, -0.5116494243639557, -0.2801781913167251, -0.3004309613612403, -0.07556511069627353, -0.39146281628711155, -0.4033401662277698, 0.38683416674421145, -0.27627061425049515, -0.36985903686769195, -0.9024072804415546, -0.47460328814609554, -0.5079150162467726, -0.4045402974103414, -0.35458196830625066, 0.1624684826797111, -0.41823695604073574, -0.47191278116518776, 100.0, -0.5504818089615208, -0.4052175396312138, -0.4636250288206107, -0.404951886891075, -0.40424242534481025, -0.29115074594485424, -0.8547757216819912, -0.6616354637330422]
    # Gini for first source = [0.41700933485324654, 0.41430692233864935, 0.4436959268797721, 0.41360228999893367, 0.5039485037059642, 0.4233504410169732, 0.4115401474295525, 0.4163730014218704, 0.5153791952172118, 0.461268819701572, 0.4967308269340974, 0.41440449916719824, 0.4187140862571708, 0.41835599517733724, 0.41493232788496287, 0.5217324124477494, 0.4886483376318146, 0.42097104330187746, 0.4188638349484116, 0.42707052316972144, 0.41963842852871414, 0.42823483249099664, 0.4161982394242443, 0.41488893256819864, 0.4192276082397819, 0.5729597574501969, 0.41348896453603, 0.4781409357531312, 0.4155438593098315, 0.41705144714993925, 0.5418633686216627, 0.4640370600630657, 0.4419343325401058, 0.4251418619695053, 0.4223582812515552, 0.45802216742162116, 0.4280224527827708, 0.451577191473874, 0.33063216809427526, 0.4137308042871799, 0.4195488444831366, 0.3059465252086053, 0.4730588838214684, 0.41920608787000635, 0.4196122153569193, 0.41645071116564875, 0.41369380439134584, 0.44369642336518256, 0.4188117805746105, 0.5066747626577905, 0.4785139780578223, 0.4297635469045746, 0.418420413580294, 0.4315761023756293, 0.5099128554951193, 0.4198764146710052, 0.4177449174799942, 0.42340769022898184, 0.5672717043020433, 0.5251160619224312, 0.42763359599036516, 0.41746626023722605, 0.41608527898843983, 0.41845148935748167, 0.4240930389994372, 0.4182398277739691, 0.42035821051256844, 0.41849923229875635, 0.41987046004985784, 0.42258123866086733, 0.42263382634683677, 0.4260644634296677, 0.4253880353385978, 0.4156192280714102, 0.4685691587265177, 0.4257765247863124, 0.5755708780765724, 0.5031676390825681, 0.4156381809772079, 0.35415422077914277, 0.5547202106229052, 0.4727876621506913, 0.4089618910696646, 0.38468971805201624, 0.41407486176783376, 0.4251157696443124, 0.45123482325943665, 0.4241800588270384, 0.5631653772911305, 0.4558890816662085, 0.4101331408178426, 0.5413699586120957, 0.4265614127768432, 0.42579429518183765, 0.48568392869373406, 0.42904917010207066, 0.46843810251032103, 0.39499430164813576, 0.31425627859156985, 0.4199504259533126]
    # M20 for first source = [-0.6915091466431722, -0.5554281967983562, -0.46394588408924203, -0.49311092340981555, 0.11789515125895017, -0.4916582814752788, -0.5419458095023758, -0.5724467695732322, -0.7076952459169502, -0.566845456369916, -0.2541311710135192, -0.5944458495034883, -0.6697514934815971, -0.5281456876722826, -0.6153649931741728, -0.6695819053120353, -0.5343614954917417, -0.6479933850401873, -0.5461116227570935, -0.6159586249985661, -0.4920361292387398, -0.5623327333389121, -0.5551565032084154, -0.5630121094686882, -0.5766256488049463, -0.6919649183610597, -0.4156062753909184, -0.6716227469185254, -0.5690652344704403, -0.5273131090801703, -0.6911054266786366, -0.3944041424369307, -0.4524705385043616, -0.5764794637977342, -0.7148101078585862, -0.5725762279583289, -0.7581815612951037, -0.4829396429055158, -0.6916860634001605, -0.6137719842120783, -0.5216168140522617, -1.1845595073473636, -0.5803560691584515, -0.5823504484889455, -0.3124642952453997, -0.5395200861868477, -0.5695224074418269, -0.5720317241484112, -0.4289365800675783, -0.7509691690884079, -0.6157525360885776, -0.5894639469799793, -0.5985664963735972, -0.6027089219266324, -0.6596883329324745, -0.6011601831647365, -0.4797891608425443, -0.48809247172105363, -0.6992002795005064, -0.48599846417820197, -0.5716007534378681, -0.5905342310433757, -0.47431745408833825, -0.6297033324995764, -0.8305255610734349, -0.5857370366673131, -0.5705508084019802, -0.5775268358244586, -0.5360611409161357, -0.5867503733783959, -0.6380472658477085, -0.8129437064359313, -0.568481934904904, -0.6951663485091076, -0.2550135555260724, -0.5327544729159365, -0.7015562159278051, -0.6108540960443293, -0.6162814104772694, -0.6928009950621311, -0.6954495100608787, -0.4898213957404206, -0.6906731267164428, -0.6915863708587371, -0.5291460008703321, -0.6120700263776957, -0.6765520444130623, -0.587966302075239, -0.7000225673518239, -0.29499076831282683, -0.5736405842063153, -0.694508700508313, -0.6560815843530248, -0.6445343859272651, 0.04429461325164261, -0.5198973667201221, -0.8328807136740619, -0.6892732683221792, -0.7168544345184584, -0.6266545794498342]

    ## for second source gs_irs45 but excluding error ones (they were regenerated)
    conc = np.array([0.6841664569043726, 2.771649759631824, 1.8192049743482746, 1.8095524739573599, 2.209011027871299, 2.220345289960115, 1.5694324433763458, 1.0661397954002567, 2.67698359496259, 1.4371639947116452, 3.1248089936473225, 1.6051510925554777, 1.4597351042682516, 1.4824078169965305, 2.2114824978152066, 0.8036525355517878, 2.0250330187048173, 3.2539147300059827, 1.7355431974365898, 1.868819857555801, 2.0251470963306764, 1.416640315212444, 1.3288637135883934, 1.505149802161684, 0.4408310272591607, 1.813337416518589, 2.5751156294963726, 3.0211062268037643, 2.7908425634887326, 2.4236625845589996, 0.8627621039652991, 1.7156523801947756, 2.3219818544833695, 1.4917613384064319, 2.471830146029438, 0.40540722047914546, 1.7655833174616422, 1.5051510678326148, 1.5289608790464062, 0.8591983100337341, 1.3296139607475586, 1.4595238911141366, 0.7877333787027643, 1.555435126439547, 1.631187940196154, 1.0726609541841505, 0.7654736879185478, 1.3706936850940972, 2.336733955595555, 4.918098914681252, 2.965444127075563, 2.5884877006882023, 1.7494204715139507, 2.7604544438795076, 2.3409984647554123, 1.3351401820224025, 0.563985431533575, 0.35352537940088313, 1.52757838635428, 1.7768934508330316, 1.146801934548575, 2.155336585030585, 0.9923232351581258, 1.8383458289284171, 0.9084606152009036, 1.6476342044237624, 2.2278731223701524, 2.772507890042883, 1.505192156807117, 1.6125666018311575, 1.7553820087928185, 2.3522474355062393, 4.36202782839737, 0.6651440388823273, 1.8738521099817984, 1.6176359268205989, 2.2162162253583753, 3.055639324089272, 2.97750299033849, 2.677544138346393, 1.2827691593072863, 2.258802205010258, 0.7851832129587278, 1.5051485607953665, 2.913652332907718, 2.9281186218642032, 2.2858770149048397, 1.4970738225977782, 1.726727161209164, 1.469273935823046, 2.059307684904382, 1.7800144012157981, 1.9432902915939225, 1.0498863915380088, 1.441373263319003, 0.5990627730699636, 2.0317575951895566, 0.6773561336534959, 1.5176503510588721, 2.171938180295496])
    asym = np.array([-0.26412401305895644, -0.3553518605157498, -0.31592943295478887, -0.4525640680856704, -0.563578441179093, -0.3136171196533975, -0.38365162126239666, -0.4593429583352259, -0.3352719980805092, -0.4129314240253162, -0.4794412664200756, -0.31799984220972477, -0.23804199190774158, -0.12353193939113039, -0.2743458449202833, -0.29645392227583073, -0.3320322540719472, -0.21198188421026185, -0.4210078775484991, -0.3054937614477966, -0.40635442635683605, -0.26369885413251354, -0.23250464573660334, -0.5193403990542029, -0.1867088749806358, -0.2518689906082243, -0.3429048407643239, -0.44205306896230073, -0.5901398112261049, -0.23212669631481236, -0.1661556829774311, -0.3836720342055978, -0.22286653365675257, -0.11850009721663589, -0.24763426519951154, -0.35481807818870237, -0.3720022433626162, -0.21139400754667922, -0.5217078305899425, -0.3524296964092212, -0.35174148176170467, -0.3625282155843904, -0.17430103043232706, -0.2791684475059393, -0.1066260330809167, -0.22159407629171032, -0.3240146670986332, -0.19383766915861336, -0.3405587313869454, -0.21744365414621095, -0.3245577356811942, -0.38400577085200555, -0.43875057244605825, -0.23815230559020772, 0.0494213415813203, -0.4742488788817984, -0.1572119296061226, -0.25315836781629425, -0.24816918470302576, -0.31430219792170083, -2.0941887341208902, -0.17845428304562438, -0.3318231567143298, -0.2930210955436104, -0.21667674665486067, -0.15969779325214628, -0.48780513648145735, -0.1940909480386674, -0.7615506972063053, -1.2982886744061959, 0.12196753896869557, -0.25766847507573454, -0.24673730475439357, -0.2741330476195321, -0.2625093494927368, -0.22789163121556827, -0.38740380771711475, -0.3748842035606513, -0.27503660629798926, -0.16054028037752716, -0.30951388771892824, -0.17626212945416134, -0.38182681675658164, -0.46527118226423025, -0.3709237717520132, -0.34223678905223043, -0.22100665967184313, -0.393173820822022, -0.32830127406066406, -0.4319942235300303, -0.39723193790211314, -0.19127261416410776, -0.23025221327359277, -0.32747119067245034, -0.17668465977401068, 0.014250904657745787, -0.08184353187105428, -0.2734152954526674, -0.8012376155864551, -0.18511490713724119])
    gini = np.array([0.4156180755835717, 0.3997975118350415, 0.42278575990378015, 0.40852629632667226, 0.35081894016594467, 0.4279755039893405, 0.4094789501196255, 0.42086613356175906, 0.4184856953214233, 0.435796828822842, 0.49380419230775824, 0.4160816766583225, 0.4164290415769248, 0.46938730949749896, 0.43318508032984554, 0.4971458650159903, 0.4146931335502107, 0.41495767489949514, 0.5537071543108982, 0.422187855346756, 0.42110288818173736, 0.41631276436351966, 0.4202737077624125, 0.4308250897511564, 0.41155521358441494, 0.4210799628692876, 0.42190203356013833, 0.45490265759898735, 0.45989284572857203, 0.4136659271937617, 0.4202546079761178, 0.38623259930480264, 0.4219988509902662, 0.47769301641117934, 0.5295489816059654, 0.422240963039649, 0.48153322750177713, 0.41430826336860127, 0.49723261131369495, 0.4301341085464227, 0.39548310416891186, 0.4171834242655077, 0.44743801419718887, 0.5247475413978067, 0.5361093279756681, 0.41820281019645733, 0.3094624491206618, 0.41779522251847295, 0.41644086769425226, 0.4600277604900386, 0.4176643672887392, 0.48820214619725144, 0.42924659943891974, 0.4159992990709674, 0.38312907824439013, 0.4134308386427213, 0.42319718146847507, 0.5315237034957477, 0.5549785565176762, 0.49238564880076596, 0.414063749062927, 0.4110264843475521, 0.4191803212296857, 0.3958860209595293, 0.41400572170669186, 0.41971487790937634, 0.40798829777930384, 0.4239636939388489, 0.41007864232799873, 0.4162726135151605, 0.4215599791636615, 0.509275249927276, 0.42114066718125887, 0.4353176618069806, 0.43430783680160534, 0.4232900889794186, 0.432113035540187, 0.4732577258278247, 0.42478721010208115, 0.5677982806184361, 0.422146816956927, 0.46522942332064904, 0.42183401941620635, 0.41240387586076055, 0.3819280463270816, 0.486013593817146, 0.4225508942103308, 0.43113689205513533, 0.4155670944528886, 0.42522562044761275, 0.41800854259894443, 0.5469511540647553, 0.42353059016709865, 0.48531886155827, 0.42040598428023834, 0.4140442000609514, 0.4349413842550235, 0.4733925040262706, 0.4352750090264443, 0.4233162375735778])
    m20 = np.array([-0.7113536102656963, -0.5553117596374736, -0.5449766189265206, -0.46389458637421743, -0.6898677318957384, -0.5603063062499726, -0.7034103981567477, -0.2918521330626047, 0.09593053822674355, -0.5426574318300572, -0.7144236134891704, -0.6150218927310338, -0.648266904793741, -0.5686998689087721, -0.3380662737158148, -0.6774524275370568, -0.8596732574853584, -0.6111939293937136, -0.6981896077279434, -0.5505812432534791, -0.6051021456103974, -0.39475070369091786, -0.5532997310741696, -0.6147598364358023, -0.4864762387684689, -0.6242706377489626, -0.6060387113943573, -0.7140082635536923, -0.6298625717405221, -0.5817016877475922, -0.5678177376901976, -0.7059934162389382, -0.6541350895168733, -0.5470518256062806, -0.6763881379100565, 1.519336220845862, -0.49553138338394964, -0.6606615748489716, -0.6858406724930757, -0.6507067966522166, -0.530006122257833, -0.5663714255448012, -0.6410636869411529, -0.5181167838056657, -0.7029778230491978, -0.6599082988669428, -0.7121398667730648, -0.5739873571465797, -0.5783754933876657, -0.7068508113244926, -0.62520905741329, -0.6973300584114298, -0.5143786266291067, -0.637497458435298, -0.688238417865544, -0.5488241840640551, -0.5636455505445132, -0.7073259655018065, -0.7038817589343853, -0.6919065620169484, -0.6483357270770533, -0.6697774633233275, -0.6858317828777681, -0.748973042016392, -0.5452866279617379, -0.4394496630307352, -0.5780979508272884, -0.6919013873174689, -0.6952380078017122, -0.6549038842489027, -0.5768859007231698, -0.4830777365611909, -0.6727828669110818, -0.5812415650543956, -0.5660669146129786, -0.48520958003121617, -0.6871281708398901, -0.621166444501068, -0.6405328767437405, -2.584175957123016, -0.6679621121656502, -0.7226937211496429, -0.4827166263305287, -0.5812079224408209, -0.5048721446467569, -0.6842682150035643, -0.4528340980761042, -0.6585442451377403, -0.5650256456251381, -0.6159184903997525, -0.607244853267156, -0.6966725662037824, -0.36875886958854803, -0.6917367419897105, -0.4152525693037835, -0.5632419266876881, -0.7526030123943591, -0.452698573021316, -0.6469517469560622, -0.5537884596293732])
    
    print(np.percentile(conc,16)-np.median(conc))
    exit()
    ## for second source GS_IRS45 and including error ones too
    # conc = np.array([1.7425965617443033, 2.4256391089642926, 0.8406873203834677, 1.852876812793498, 2.1280945575046784, 1.71221601500097, 1.159779520182173, 5.229969983978497, 1.304842095354814, 1.7756520515920382, 1.4635871926312014, 3.1432044295267128, 1.035421567038611, 1.0034399867010615, 2.6497865039563036, 1.505150957788011, 1.5778758371135688, -99.0, 1.7618274011551978, 1.4152393284124258, 1.479606835196674, 1.667706168340637, 2.298071603222974, 0.8072106082442465, 2.3806861157520505, 1.8048443517837265, 0.29571864285904714, 0.6931934469535348, 2.1686036774108564, 1.5167101911210665, 1.9689423846528322, 0.542721110836813, 0.0, 1.24233268544634, 1.6181712928724075, 1.46304836927114, -99.0, 1.2638363456178734, 2.815753708252005, 1.7548293781719169, -99.0, 1.3005548492027141, 2.0860990511159296, 2.0843988303077805, 3.7507679206063904, 2.291053765132076, 1.7773148757164647, 2.3794311681596376, 1.787280590721678, 2.1251293428128792, 1.566724972926819, 2.8383069864230395, 1.6263310023684574, 2.5512506565993673, 1.1905447049595983, 0.8893143080824001, 3.607326775459581, 2.5025083033850315, 3.133737171838633, 1.607545501011015, 2.1000120277832974, 1.759600315324994, 1.948202579874216, 0.9125073691850749, 0.31639263187768113, 1.4604027693836519, 1.8068544603692718, 2.183071201031007, 2.4498877013243443, 2.0519034896122172, 3.3512831688739193, 2.6990561781696805, 2.0394114322439516, 1.4557662372950149, 1.9899353090882337, 1.4589930150468686, 1.6623732044763573, 1.9204134255357963, 1.5624310052129644, 2.730670106105148, -99.0, 2.561784147145918, 1.7205002919576868, 1.5785182412587146, 1.8356084350304527, 2.4191297168556862, 1.834020310546884, 1.9560457574056782, 1.3237368468500765, 1.5935708249671432, 1.2743264534546022, 1.3153650923481743, 1.2389069074330585, 0.0, 1.5051499764910306, 1.6274851902670213, 1.684483527262182, 2.0155606940769695, 1.469618749841565, 4.07379133250894])
    # asym = np.array([-0.1601400018194404, -0.3793838471222359, -0.22931328908412088, -0.37653050823800055, -0.263046185322643, -0.16534125056687768, -0.1310679078092471, -0.3587336875226504, -0.27782085621075886, -0.17431222641118918, -0.4143374122735252, -0.23050813149896227, -0.36413536815587116, -0.23151485582408476, -0.17101130198694509, -0.3834997475191207, -0.19740474486564485, -99.0, -0.2993868013935948, -0.2634313551304939, -0.7535200190812744, -0.3671766558805829, -0.3710144672803423, -0.1313263005066896, -0.4477576674467558, -0.42737356659054887, -0.258317824928531, -0.33123078336513934, -0.2746971598149606, -0.4322130112570039, -0.2817182819430103, -0.09518068756396093, 100.0, -0.3025329367167217, -0.38374974852925026, -0.30636151417573565, -99.0, -0.10915776326509169, -0.58308706552319, -0.5232483222021305, -99.0, -0.26125678245445505, -0.300411873734684, -0.3772061361183323, -0.293465182545479, -0.17179178611133472, -0.20062800969164638, -0.300749487859739, -0.7587362448451949, -0.40554609990984813, -0.6370817596251794, -0.22104505183277928, -0.2059250492970027, -0.3954183959325909, -0.20938520668041408, -0.3507164664464594, -0.27243326884792, -0.25373376276173604, -0.2688404569256348, -0.20301825602327903, -0.4672292829746799, -0.3696477129689131, -0.4368353596205505, -0.22926590114328801, -0.3934372110019975, -0.29097303094083676, -0.2921774728036024, -0.9225030560015309, -0.1154037855507471, -0.18779918949841967, -0.3239915353539913, -0.33440005678320167, -0.3399437005467907, -0.17785334151378698, -0.2623633699884605, -0.29285602502279307, -0.34415535108195394, 0.020497758100825273, -0.3999886470571082, -0.7114837743085682, -99.0, -0.6927843856178675, -0.337340719032557, -0.6323126620510945, -0.2694692859421593, -0.3797644707629147, -0.14378557054052643, -0.28887681724167097, -0.15393650695418964, -0.16858663459243461, -0.3283380389381548, -0.27929783508933514, -0.2051047647049347, 100.0, -0.39294655165053616, -0.1251650781527155, -0.5599771179713422, -0.30521714798856164, -0.24622543210796902, -0.3134401558491968])
    # gini = np.array([0.4319853589556525, 0.4190604476599241, 0.4311627657664521, 0.4338482480198404, 0.4150845712572121, 0.48274202042311193, 0.4929508134556021, 0.4826853666436239, 0.46388245154210955, 0.43380563362187197, 0.5125846069494966, 0.416824915788275, 0.5096998599456015, 0.41509566768895295, 0.41459494297004196, 0.416293745156061, 0.42075309060863536, -99.0, 0.5083423134190914, 0.41430772975850416, 0.41833553838721216, 0.41087553992410836, 0.4673341751582806, 0.4813875575219897, 0.44493479924233026, 0.43508267025423664, 0.41881532790010356, 0.29147317356737307, 0.41863733122600777, 0.42934456987646125, 0.41671351763371023, 0.41913576201115393, 0.5126498132053313, 0.41931743674313887, 0.5290703690514009, 0.4148813537886092, -99.0, 0.41447696982993343, 0.40032241952927666, 0.4240892006935013, -99.0, 0.5216870394891054, 0.41345124689160795, -99.0, 0.4243135283454956, 0.27362188947222293, 0.4127512367781638, 0.49852077210146145, 0.4313871815381681, 0.4199021360740269, 0.4126043413124923, 0.4740174227052421, 0.42008864184895295, 0.41834654918799014, 0.4676513433577353, 0.42334289126207525, 0.42036490684235966, 0.418698583251313, 0.41872975414182156, 0.44455017407554986, 0.40774544552079567, 0.3672824410129, 0.3462895011013332, 0.42108229943895864, 0.4541025655446262, 0.48616870638456755, 0.43005210930814686, 0.42027701919731486, 0.5117863897514853, 0.4210915298752008, 0.41947769476282765, 0.41840448426082016, 0.41516512789924237, 0.40677318089600173, 0.41496333533878066, 0.4204416955987918, 0.4738317510682166, 0.39270369600548005, 0.5425882023263922, 0.4191248543020174, -99.0, 0.42046114581108535, 0.42146772998238125, 0.49140800573807286, 0.38154992817398936, 0.4161856012804884, 0.46999572651980576, 0.4288389031729174, 0.414800919767121, 0.47891855240073955, 0.4873541400160723, 0.4198691187107796, 0.41432005117630094, 0.5158368276804011, 0.3517694994988194, 0.5192470682657802, 0.4383115957393175, 0.4239315564392656, 0.4664886119154924, 0.4610190042578979])
    # m20 = np.array([-0.7108225467611131, -0.6009518685447532, -0.25794102589041135, -0.7935063035818873, -0.5555833997392836, -0.6522205483675111, -99.0, -0.538582621391467, -0.7511586151900022, -0.5933221494645914, -0.36690391024679647, -0.5549031933287286, -0.7003803105240268, -0.669163104190716, -0.4618539181077544, -0.5311427189095426, -0.46602385556085746, -99.0, -0.6419892939726705, -0.6086941232286285, -0.5984948194011807, -0.4797407074810726, -0.6105168154520545, -0.6879473527684363, -0.5973296083833275, -0.4954932479203869, -0.6684572709222337, -0.4007624556699298, -0.6430183214284233, -0.7235986680594367, -0.8371207725751573, -0.6788539493777527, -0.6859392657851637, -0.5851432231089736, -0.689460766214362, -0.4533455403336577, -99.0, -0.5207813725150895, -0.6816710629387877, -0.4051546257797494, -99.0, -0.6412098591823577, -0.5875781105312865, -99.0, -99.0, -0.7506837366197965, -0.5481347952209416, -0.8291433997317609, -0.6152901448968736, -0.4889768104468751, -0.5242567632853928, -0.6881459925429182, -0.5840155126365499, -0.6127127810093871, -0.6611209810772994, -0.6617897219425913, -0.5527079868026753, -0.6128959331980379, -0.6177990599507477, -0.49085914797504926, -0.6811184501628057, -0.6344323689548248, -1.3649266169483176, -0.6120685465997998, 0.32606121703135504, -0.6740476759438656, -0.8108457577560955, -0.6560678969811102, -0.7079666626982364, -0.6050756191967727, -0.5852214208400305, -0.4914627196059726, -0.5869142523304944, -0.6005675705412846, 0.004614293574427743, -0.606857920290297, -0.4046519629233732, -0.6920196954558716, -0.7056693598968123, -0.5516028647373947, -99.0, -0.49620173993221767, -0.7184531194087468, -0.5128642410732607, -0.6938897845238466, -0.5001565395317642, -0.6582130981885347, -0.037059034951535175, -0.48753529191740574, -0.404106951061712, -0.7032795544601871, -0.40145592910720884, -0.7351477755317893, -0.6982759027540986, -0.7062374000515796, -0.7017506621985937, -0.5239182031843622, -0.5454553213609009, -0.5192059781869837, -0.5703250397523292])
    
    ## for first source gs_irs12 and including error ones too
    # conc = np.array([3.433982145996324, 0.0, 0.0, 2.135488498034432, 2.1364234821626944, 1.8350153492362886, 1.6648924783929173, 3.2642312598353524, 0.0, 2.867822123340811, 1.68804873589931, 2.9434437132765168, 1.0185669059846334, 1.000783408202949, 1.7639501441697734, 2.857098243903114, 1.8752388437390555, 1.591077742245046, 3.1422985679624253, 0.9936209571204959, 2.7260106374780064, 3.2411750862358932, 0.7408190248112008, -99.0, 4.088763433571802, 2.3931185464198457, 0.7651054338987175, 1.3709267963432337, -99.0, 0.0, 2.9720390094817817, 0.7865898293187905, 1.1904150237353333, 3.967126164673091, 1.7325332380529923, 1.8254731710483771, 2.9972951433971344, 1.707424487972977, 2.667542564042469, 1.153479013319993, 2.84147731065455, 3.535197707543518, 1.6841884576472597, 1.8399851594954053, 1.244540721977654, 0.23359034810374799, 1.7659706845191825, 1.6550999727660842, 1.7389280349543075, 2.9344066947973184, 1.6426614477755703, 3.9394218688568956, 0.8210317469550876, 1.938598544360886, 2.07799730192548, 2.9161973558493823, 1.2646390723708165, 0.6219864271063831, 1.1634214411003132, 1.12088750364417, 1.065047485107617, 2.801914955442961, 1.6238565826919458, 0.9350957478728362, 1.1196186999590931, 1.372662944507671, 1.709441334514164, 1.4766085012751917, 0.8783721765344578, 2.966862199019328, 0.6702822868480047, 1.280456176713253, 1.000023993453164, 1.4465353615216203, 3.9688661360284705, 2.841327166105392, 1.2452850216083098, 2.441556577697637, 1.6291097663917657, 2.755377094099136, 1.051677593399467, 2.852885011876072, 1.524085998616518, 2.736026316771457, 1.9209218519342908, 0.9763685812110517, 3.0579223989067135, 1.5345043867556851, 2.413646543312998, 1.8783863482717784, 2.933920316261859, 0.8629715624475276, 1.3648347832742838, 0.8206758030033459, 1.3581961894149681, 0.22726447880776282, 1.208226529638531, 0.7028928851184253, 0.6173198337995396, 0.6913575489797251])
    # asym = np.array([-0.579856737092189, -0.7565197719409, 100.0, -0.6589202497583144, -0.588748436800838, -0.7545713207363275, -0.9976400097561896, -0.49509676139491554, 100.0, -0.530589636484486, -0.47217113244950304, -0.5478498462409291, -0.6380175336375737, -0.6442234587488647, -0.7754706381199999, -0.6745514747633248, -0.544543900133235, -0.4617929920766424, -0.47984440477343415, -0.4338307975454217, -0.3085195119407209, -0.526645341998763, -0.38031400304965773, -99.0, -0.9359879719363914, -0.5440123606920932, -0.6964368825060491, -0.36082166162830814, -99.0, -0.6358986018680314, -0.14327809339403164, -0.4012031215085022, -0.6773868551822892, -0.6434872016395154, -1.0823566859982379, -0.40390663421196005, -0.3035125499642112, -0.4667554365183227, -0.27735685572074026, -0.3773812521129658, -0.5768474227405321, -0.5680646687008286, -0.6302221128318849, 0.05946395289576541, -0.3933242126246976, -0.948875345903831, -0.5073413335511718, -0.5177545999832195, -1.0287844930155223, -0.7011618591694182, -0.20045331468141184, -0.31605917995071253, -0.8527746544325691, -0.4267698416756679, -0.3224516965915787, -0.5538491141681696, -0.4538191100726038, -0.9102585797262711, -0.39032078812157933, -0.4090065201751203, -0.4934702190628658, -0.4476423442546901, -0.9877099453449889, -0.7947126475990126, -0.6778083996225374, -0.40664109061555376, -0.6487246490939106, -0.581631256576247, -0.44865284541612727, -0.470859766234475, -0.46534679406442897, -0.4128920403070581, -0.4615525123364108, -0.4193824338619551, -0.5024197463273531, -0.400092170017443, -1.1721287318962013, -0.6886207339656898, -0.3955955408650267, -0.3324159815978434, -0.325706463765071, -0.43827997166102084, -0.5483205137058011, -0.5753129816814054, -0.3161595945359515, -0.5640848990710716, -0.0479654998583268, -0.5900606518731886, -0.5087813350755261, -0.1600167658845579, -0.5475909352338971, -0.6268802488623134, -0.21210428856860633, -0.47636044154095597, -1.11682539178164, -0.3921873323065164, -0.45770445464646997, -0.5431551651313994, -0.27675605237235945, -0.3344020031311419])
    # gini = np.array([0.4245272940438047, 0.44597881740665385, 0.6352719582523028, 0.6495436766546487, 0.47946276566820634, 0.43240630969572036, 0.5449404265020593, 0.5497086596054779, 0.6357633496582956, 0.48401856778770724, 0.48329184665523933, 0.5803354752783858, 0.5735818195357534, 0.43882971066962784, 0.558458786542014, 0.6174401858111169, 0.4804059902824454, 0.4929826771520554, 0.5672111048054953, 0.462548571593222, 0.5582973355175855, 0.5090592273873323, 0.4479444967661722, -99.0, 0.5339009176197772, 0.5582354699579246, 0.4736898861350359, 0.43298470107198717, -99.0, 0.4279247740091807, 0.4401854377113145, 0.4305132983190205, 0.4957242827167998, 0.4557362900921032, 0.41825059261904635, 0.44027420095049524, 0.5739698099166011, 0.44322656627571755, 0.6305105217079832, 0.4610659530197838, 0.44062995263486926, 0.46777224787512967, 0.4600955870826492, 0.35290520839560696, 0.4438239758469482, 0.4701729353720021, 0.4407322979491228, 0.4796344058666678, 0.4331518553498912, 0.5363701545634657, 0.39767816203662815, 0.6592480535294362, 0.5395885797468188, 0.5628145652812175, 0.40402000482511274, 0.449496252645953, 0.5686486670121245, 0.44788488160439355, 0.4257397606726324, 0.47818920366763823, 0.5459445647283787, 0.5177197674825031, 0.5242924859558534, 0.4423322506431304, 0.6002612864954412, 0.46715488196760846, 0.3242047391185894, 0.4533868713875001, 0.3373637787873169, 0.49095823805825123, 0.4392712049181762, 0.47413818833851046, 0.5430549900614725, 0.4642888150107608, 0.6190759123986898, 0.45565636615232685, 0.5359593651947803, 0.6161568389869098, 0.47030980281939166, 0.44219616767795406, 0.6504632826373234, 0.511601739488806, 0.6498100846385932, 0.5565128686322063, 0.44102136351603916, 0.6165623197296639, 0.4540076013231611, 0.5083234020561413, -99.0, 0.6565292127846594, 0.49703417920662685, 0.44122130880998695, 0.4283705477291366, 0.42588092199597183, 0.46941086418334144, 0.44614295861165093, 0.5038646060598783, 0.4706700125032189, 0.4310409165992756, 0.4452111869929273])
    # m20 = np.array([-1.6032438035480983, -1.7022012682469194, -0.7398066259343489, -0.7318460562811874, -0.8444879971700647, -1.2781262234690876, -1.1391981120234975, -1.5319125401592069, -0.744057203246009, -1.9515490462504648, -0.9276003591171025, -1.336168717106628, -1.6753679499068959, -0.8118339601797524, -1.1464712332593907, -0.7442029637648173, -1.4844429176766027, -1.5683434166944068, -1.506731626136148, -1.712354247745691, -0.833554392418617, -1.5484298830954315, -1.472488592100295, -99.0, -1.1282623058555514, -99.0, -0.7521412839763575, -0.7335516117141415, -99.0, -1.7810909845913312, -1.0102256154280764, -0.7575448531624867, -0.9240790131248974, -0.8585489178593766, -0.7390823517080783, -0.7726143965281551, -0.8368746027119082, -1.6292184518061892, -0.7488180446003732, -0.7629739333105796, -0.7287712211596221, -1.5637317784995581, -0.7703333054476706, -0.7337470461733467, -1.4026119679063702, -1.8387651884740002, -0.9518266603510774, -1.0132379638842028, -0.8867221673428567, -1.7349515339844954, -1.2142594136298452, -0.7472548531070843, -1.1743223344084992, -1.5249001676505922, -0.7567801980729303, -0.7233440941665019, -1.6562499478506363, -1.0794306106122626, -0.9287270042736732, -0.8857621284938159, -0.8962193748168671, -1.0536763524092305, -1.671698111624322, -2.192454192376976, -0.7985974856523695, -1.0266776498102566, -0.7937145453848069, -1.2010279132102009, -0.7371392773566222, -1.354470817698223, -0.9851069670651293, -1.9852685632920155, -1.2791246741966154, -1.48152040064554, -0.7865541735302871, -1.2158033045932002, -1.5038723899512492, -0.792180642759663, -1.2660837491818913, -1.0651490509815944, -0.7412629521511531, -1.7607689692319337, -0.7475992786213924, -1.6347709544787015, -0.7526022031084749, -0.754738542119239, -0.9489823521459767, -0.8981166983929195, -99.0, -0.7539002017623373, -0.8658701830435247, -1.7190961678188652, -0.7867937723468862, -1.2474460654135475, -0.8162037072786045, -0.9932128723707554, -1.5902350592775283, -0.836001912845123, -0.9578679351250263, -1.7169333634842892])

    histbins = 7
    # bins=np.histogram(np.hstack((filtered_group.get_group('f150w')['ID'],filtered_group.get_group('f200w')['ID'],filtered_group.get_group('f277w')['ID'],filtered_group.get_group('f356w')['ID'],filtered_group.get_group('f444w')['ID'])), bins=histbins)[1] #get the bin edges

    # ax[0,0].hist(filtered_group.get_group('f150w')['Concentration (C) Error (84%)']-filtered_group.get_group('f150w')['Concentration (C) Error (16%)'],bins=5,histtype='step')
    # ax[0,0].hist(filtered_group.get_group('f200w')['Concentration (C) Error (84%)']-filtered_group.get_group('f200w')['Concentration (C) Error (16%)'],bins=5,histtype='step')
    # ax[0,0].hist(filtered_group.get_group('f277w')['Concentration (C) Error (84%)']-filtered_group.get_group('f277w')['Concentration (C) Error (16%)'],bins=15,histtype='step')
    # ax[0,0].hist(filtered_group.get_group('f356w')['Concentration (C) Error (84%)']-filtered_group.get_group('f356w')['Concentration (C) Error (16%)'],bins=5,histtype='step')
    # ax[0,0].hist(filtered_group.get_group('f444w')['Concentration (C) Error (84%)']-filtered_group.get_group('f444w')['Concentration (C) Error (16%)'],bins=5,histtype='step')
    for i in filtered_group:
        ### for now it's just the histogram of a few sources, will eventually do the full thing but it'll take a while to run all them and save them
        ### idk if I can save an np array to a table like that        
        if(i[0]!='f277w'):
            bins = 5
        else:
            bins = 20
        # ax[0,0].hist(conc-np.median(conc),5,alpha=0.5,label=i[0])
        ax[0,0].hist(i[1]['Concentration (C) Error (84%)']-i[1]['Concentration (C) Error (16%)'],bins=bins,label=i[0],histtype='step')
        # # # plt.xlim([-10,10])
        ax[1,0].hist(i[1]['Asymmetry (A) Error (84%)']-i[1]['Asymmetry (A) Error (16%)'],bins=bins,label=i[0],histtype='step')
        # ax[1,0].hist(asym-np.median(asym),5,alpha=0.5,label=i[0])
        ax[1,1].hist(i[1]['Gini Error (84%)']-i[1]['Gini Error (16%)'],bins=bins,label=i[0],histtype='step')
        # ax[1,1].hist(gini-np.median(gini),5,alpha=0.5,label=i[0])
        ax[0,1].hist(i[1]['M20 Error (84%)']-i[1]['M20 Error (16%)'],bins=bins,label=i[0],histtype='step')
        # ax[0,1].hist(m20-np.median(m20),5,alpha=0.5,label=i[0])

    ax[0,0].set_xlabel('Concentration')
    ax[0,0].legend()
    ax[1,0].set_xlabel('Asymmetry')
    ax[1,1].set_xlabel('Gini')
    ax[0,1].set_xlabel('M20')
    plt.show()

error_dist_hist()
exit()


    


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

def A_S_scatter():
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

def A_C_scatter():
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
            # plt.legend()
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

A_agn_scatter()
exit()

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

A_agn_histogram()
exit()

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
        plt.errorbar(x=(np.max(i[xsearch])+np.min(i[xsearch]))/2,y=np.median(i[ysearch]),yerr=np.std(i[ysearch]),c='gray',fmt='.',alpha=0.5)

    # plot the scattered values
    ind = 0
    # for g in [agn_class,comp_class,sf_class]:

    for c in main_colors:
        try:
            # y_err = np.array([combined_class.get_group(c)[f'{ysearch} Error (16%)'],combined_class.get_group(c)[f'{ysearch} Error (84%)']])
            # plt.errorbar(x=combined_class.get_group(c)[xsearch],y=combined_class.get_group(c)[ysearch],yerr=y_err,fmt=',')
            plt.scatter(combined_class.get_group(c)[xsearch],combined_class.get_group(c)[ysearch],40,marker=main_colors[c][1],facecolor="none", 
                        edgecolor=main_colors[c][0],lw=1.5,label=(f'{c}'))
            # plt.legend()
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

def S_agn_scatter():
    ### not using as of 7/29 ------

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
        plt.errorbar(x=(np.max(i[xsearch])+np.min(i[xsearch]))/2,y=np.median(i[ysearch]),yerr=np.std(i[ysearch]),c='gray',fmt='.',alpha=0.5)


    # plot the scattered values
    ind = 0
    # for g in [agn_class,comp_class,sf_class]:
    for c in main_colors:
        try:
            plt.scatter(combined_class.get_group(c)[xsearch],combined_class.get_group(c)[ysearch],40,marker=main_colors[c][1],facecolor="none", 
                        edgecolor=main_colors[c][0],lw=1.5,label=(f'{c}'))
            # plt.legend()
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
    bins=np.histogram(np.hstack((agn_dom[xsearch],comp_dom[xsearch],sf_dom[xsearch])), bins=histbins)[1] #get the bin edges
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
        
        plt.hist(g[xsearch],bins=bins,label=lab[ind],color=color_agn[ind],fill=filling,histtype='step',linewidth=2)
            # plt.legend()
        plt.axvline(np.median(g[xsearch]),color=color_agn[ind],linestyle='--',linewidth=2)
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

gini_agn_histogram()
exit()


def blank_agn_plot():
    # template for making the cutout plot of type and agn
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
    color_agn = ['#1f77b4','#ff7f0e','#2ca02c']
    
    
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
    plt.scatter(sf_dom[xsearch],agn_dom[ysearch],label=lab[0],color=color_agn[0])
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
                plt.scatter(g.get_group(c)[xsearch],g.get_group(c)[ysearch],40,marker=main_colors[c][1],facecolor="none", edgecolor=color_agn[b],lw=1.5,label=(f'{lab[b]}-{c}'))
                # plt.legend()
            except:
                continue
        b+=1
    # plt.legend()
    # plt.show()

    highlight_mergers(plt,x=xsearch,y=ysearch)
    highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)

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

gini_m20_scatter()
exit()

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

def gini_A_scatter():
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
        plt.errorbar(x=(np.max(i[xsearch])+np.min(i[xsearch]))/2,y=np.median(i[ysearch]),xerr=np.std(i[xsearch]),yerr=np.std(i[ysearch]),c=color_agn[b],fmt='.',alpha=0.5)
        b+=1


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
    highlight_clumps_and_arms(plt,x=xsearch,y=ysearch)

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

m20_A_scatter()
exit()


#### older graphs that haven't been update (7/15)
# type_z_hist()

# testing_agn_frac_z()

# agn_frac_merger_scatter()

# agn_frac_z()
######################################################
#### used graphs as of (7/15)
# agn_frac_hist()
# agn_frac_merger_hist()
# compare_subsample()
## non parametric measurements

# A_S_scatter() # not using this one as of 7/29
# A_C_scatter()
# A_agn_scatter()
# C_agn_scatter()
# S_agn_scatter() # not using this one as of 7/29
# gini_m20_scatter()




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