import matplotlib.pyplot as plt
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




##### get the merger variants for each classification
##### for now we're including the 'maybe' flag, but in the future we should exclude it (or just exclude if from the sheet)
# v = 0
mergers = np.empty(len(morph_types),dtype=object)
mergers = table.query('`Merger (flag threshold)` == "yes"')
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
print(f'Total # of sources in original table: {table.shape[0]}')
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


base_wv_rest = table['Rest WV (um)'].mean()
err_wv_rest = table['Rest WV (um)'].std()

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
    plt.title('$\lambda_{rest}$ Distribution ($\mu m$)')
    plt.xlabel('$\lambda_{rest} (\mu m)$')
    plt.legend([f'$\lambda_{{rest}}={base_wv_rest:.2}\pm{err_wv_rest:.2}$'])
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

    agns = table.groupby(['AGN (bin)','Classification']).size().unstack()
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
    sf_sample_table = table.query('`AGN(%)`<20')
    comp_sample_table = table.query('`AGN(%)`>=20 & `AGN(%)`<80')
    pure_sample_table = table.query('`AGN(%)`>=80')
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
    table.insert(2,"Sub-sample",df['Sub-sample'].values)
    count_data = table.groupby(['Classification', 'Sub-sample']).size().unstack(fill_value=0)
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

### import the statmorph measurements tsv files ***
### only importing one filter for now just for testing
# nir_filter = 'f356w'
# nir_filters = ['f277w','f356w','f444w']
# working with grizli measurements for now
stat_measures = f"research/statmorph_output/grizli/grizli-rest-filters.tsv"
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


# combined_class = combined.groupby('Classification')
combined_class = combined.groupby('Main Class')
all_mergers = combined.query('`Merger (flag threshold)`=="yes"')

def highlight_mergers(plot,x='AGN(%)',y='AGN(%)'):
    plot.scatter(all_mergers[x],all_mergers[y],marker=11,label='Mergers',color='black')
    return plot




##############################################
####### taken from Conselice 2014   ##########
####### this is for the optical R-band   #####
# conselice_measurements = pd.read_csv(f'{working_directory}conselice-measurements.csv',sep=',')
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


def A_S_scatter():
    ### not using as of 7/29
    
    ### scatter the asymmetry vs smoothness (clumpiness) measurements

    fig = plt.figure(figsize=(10,6))

    ## overplot the Conselice values in gray to show the different ranges
    # conselice_measurements.plot(x='Clumpiness (R)',y='Asymmetry (R)',kind='scatter',xerr='Clumpiness (+/-)',yerr='Asymmetry (+/-)',c='gray')
    for g in conselice_values:
        clumpi = conselice_values[g][4]
        clumpi_err = conselice_values[g][5]
        asymm = conselice_values[g][2]
        asymm_err = conselice_values[g][3]
        plt.errorbar(clumpi,asymm,xerr=clumpi_err,yerr=asymm_err,fmt='.',alpha=0.5,color='darkgray')
        plt.annotate(g,(clumpi,asymm))
    # plt.scatter(agn_dom['Smoothness (S)'],agn_dom['Asymmetry (A)'],label='AGN')
    # plt.scatter(comp_dom['Smoothness (S)'],comp_dom['Asymmetry (A)'],label='Composite')
    # plt.scatter(sf_dom['Smoothness (S)'],sf_dom['Asymmetry (A)'],label='SF')

    lab = ['AGN','Composite','SF']
    color_agn = ['#1f77b4','#ff7f0e','#2ca02c']

    b=0
    for g in [agn_class,comp_class,sf_class]:
        for c in bin_colors:
            try:
                plt.scatter(g.get_group(c)['Smoothness (S)'],g.get_group(c)['Asymmetry (A)'],marker=bin_colors[c][1],color=color_agn[b],label=(f'{lab[b]}-{c}'))
                # plt.legend()
            except:
                continue
            # plt.scatter(g.get_group('Mostly Spheroid')['Smoothness (S)'],g.get_group('Mostly Spheroid')['Asymmetry (A)'],marker=bin_colors['Mostly Spheroid'][1],color=color_agn[b],label=lab[b])
        b+=1
    # plt.legend()
    # plt.show()

    highlight_mergers(plt,x='Smoothness (S)',y='Asymmetry (A)')
    # plt.scatter(all_mergers['Smoothness (S)'],all_mergers['Asymmetry (A)'],marker=11,label='Mergers',color='black')
    plt.xlabel('Clumpiness (S)')
    plt.ylabel('Asymmetry (A)')
    plt.xlim([-.01,.25])
    plt.legend(loc='best',prop={"size":8})
    plt.tight_layout(pad=2)

    plt.title(f'{nir_filter.upper()} A vs S')
    plt.gca().invert_yaxis()
    # plt.savefig(f'{working_directory}/Asymmetry_Clumpiness/A_S_scatter_with_conselice_full_{nir_filter}.png',dpi=300)
    plt.show()
    plt.close()

def A_C_scatter():
    ### scatter the asymmetry vs smoothness (clumpiness) measurements
    xlabel = 'Concentration (C)'
    ylabel = 'Asymmetry (A)'
    savename = f'/Asymmetry_Concentration/A_C_conselice_full.png'
    xsearch=''
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))

    ## overplot the Conselice values in gray to show the different ranges
    for g in conselice_values:
        conc = conselice_values[g][0]
        conc_err = conselice_values[g][1]
        asymm = conselice_values[g][2]
        asymm_err = conselice_values[g][3]
        if(g=='Ellipticals'): # label it 
            plt.errorbar(conc,asymm,xerr=conc_err,yerr=asymm_err,fmt='.',alpha=0.5,color='darkgray',label='Conselice 2014')
        else: # dont label
            plt.errorbar(conc,asymm,xerr=conc_err,yerr=asymm_err,fmt='.',alpha=0.5,color='darkgray')
        plt.annotate(g,(conc,asymm),alpha=.8)
    
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
                # print(f'errrrrrrrr ----- on classification "{c}" for agn class: {b}')
                continue
            # plt.scatter(g.get_group('Mostly Spheroid')['Smoothness (S)'],g.get_group('Mostly Spheroid')['Asymmetry (A)'],marker=bin_colors['Mostly Spheroid'][1],color=color_agn[b],label=lab[b])
        b+=1
    # plt.legend()
    # plt.show()
    highlight_mergers(plt,x=xsearch,y=ysearch)
    # plt.scatter(all_mergers['Smoothness (S)'],all_mergers['Asymmetry (A)'],marker=11,label='Mergers',color='black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.xlim([-.01,.25])
    plt.legend(loc='lower right',prop={"size":8})
    plt.tight_layout(pad=2)

    plt.title(f'{ylabel} vs {xlabel}')
    plt.gca().invert_yaxis()
    plt.savefig(f'{working_directory}/{savename}',dpi=300)
    plt.show()
    plt.close()
    # exit()

# A_C_scatter()

def A_agn_scatter():
    ### scatter the given xlabel vs ylabel & include highlighting of different visual classifications
    ### give savename (include file extension (.png))
    ### xsearch & ysearch are used for the table searching (if it's different than what you're labeling it as)
    ### also highlight the agn regions (agn, composite, sf)
    xlabel = 'AGN (%)'
    ylabel = 'Asymmetry (A)'
    savename = f'/CAS_vs_AGN/A_AGN_full.png'
    xsearch='AGN(%)'
    ysearch=''

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(8,6))


    # plot the conselice asymmetry values
    # plt.plot(x, y, 'k-')
    # plt.fill_between(x, y-error, y+error)

    # plot the scattered values
    ind = 0
    # for g in [agn_class,comp_class,sf_class]:
    for c in main_colors:
        try:
            plt.scatter(combined_class.get_group(c)['AGN(%)'],combined_class.get_group(c)['Asymmetry (A)'],marker=main_colors[c][1],label=(f'{c}'))
            # plt.legend()
        except:
            print(f'err with a_agn_scatter() for main_color: {c}')
            continue
    ind+=1

    # highlight the mergers (if any)
    highlight_mergers(plt,x='AGN(%)',y='Asymmetry (A)')


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
    # plt.xlim((-.025,.125))
    plt.title(f'{ylabel} to {xlabel} Comparison')
    plt.savefig(f'{working_directory}/{savename}')
    plt.gca().invert_yaxis()
    plt.show()
    plt.close()

def C_agn_scatter():
    ### scatter the given xlabel vs ylabel & include highlighting of different visual classifications
    ### give savename (include file extension (.png))
    ### xsearch & ysearch are used for the table searching (if it's different than what you're labeling it as)
    ### also highlight the agn regions (agn, composite, sf)
    xlabel = 'AGN (%)'
    ylabel = 'Concentration (C)'
    savename = f'/CAS_vs_AGN/C_AGN_full.png'
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

    # plot the scattered values
    ind = 0
    # for g in [agn_class,comp_class,sf_class]:
    for c in main_colors:
        try:
            plt.scatter(combined_class.get_group(c)[xsearch],combined_class.get_group(c)[ysearch],marker=main_colors[c][1],label=(f'{c}'))
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
    plt.tight_layout(pad=2)
    plt.title(f'{ylabel} to {xlabel} Comparison')
    plt.savefig(f'{working_directory}/{savename}')
    # plt.gca().invert_yaxis()
    plt.show()
    plt.close()

def S_agn_scatter():
    ### not using as of 7/29 ------

    ### scatter the given xlabel vs ylabel & include highlighting of different visual classifications
    ### give savename (include file extension (.png))
    ### xsearch & ysearch are used for the table searching (if it's different than what you're labeling it as)
    ### also highlight the agn regions (agn, composite, sf)
    xlabel = 'AGN (%)'
    ylabel = 'Clumpiness (S)'
    savename = f'/CAS_vs_AGN/{nir_filter}_S_AGN_conselice_full.png'
    xsearch='AGN(%)'
    ysearch='Smoothness (S)'

    xsearch=xlabel if xsearch=='' else xsearch
    ysearch=ylabel if ysearch=='' else ysearch
    
    fig = plt.figure(figsize=(10,6))


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
    plt.annotate('Merger',(-.75,0.65),color='darkgray')

    ## show E/S0/Sa & Sb/Sc/Ir separator line (Lotz et al. 2008)
    x = [-1.7,-3.0]
    y = [0.568,0.38]
    plt.plot(x, y,color='darkgray')
    plt.annotate('E/S0/Sa',(-2.5,0.55),color='darkgray')
    plt.annotate('Sb/Sc/Ir',(-1.2,0.42),color='darkgray')

    # plt.scatter(agn_dom[xsearch],agn_dom[ysearch],label=lab[0],color=color_agn[0])
    # plt.scatter(comp_dom[xsearch],comp_dom[ysearch],label=lab[1],color=color_agn[1])
    # plt.scatter(sf_dom[xsearch],sf_dom[ysearch],label=lab[2],color=color_agn[2])


    b=0
    for g in [agn_class,comp_class,sf_class]:
        for c in main_colors:
            try:
                plt.scatter(g.get_group(c)[xsearch],g.get_group(c)[ysearch],marker=main_colors[c][1],color=color_agn[b],label=(f'{lab[b]}-{c}'))
                # plt.legend()
            except:
                continue
        b+=1
    # plt.legend()
    # plt.show()

    highlight_mergers(plt,x=xsearch,y=ysearch)

    # plt.scatter(all_mergers['Smoothness (S)'],all_mergers['Asymmetry (A)'],marker=11,label='Mergers',color='black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.xlim([-.01,.25])
    plt.legend(loc='best',prop={"size":9})
    plt.tight_layout(pad=2)

    plt.title(f'{ylabel} vs {xlabel}')
    plt.gca().invert_xaxis()
    plt.savefig(f'{working_directory}/{savename}',dpi=300)
    plt.show()
    plt.close()


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
A_agn_scatter()
C_agn_scatter()
# S_agn_scatter() # not using this one as of 7/29
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