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
ti_m = os.path.getmtime(file)
m_ti = time.ctime(ti_m)

def disclaimer():
    print(f"[!!!] Reference file was last updated at: {m_ti}")

disclaimer()


table = pd.read_csv(file,sep='\t')
table.sort_values(by='Properties',ascending=True,inplace=True)

def remove_excluded(tab,name='Properties'):
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
    disk_only = tab.query('`Mostly Disk (4.3 Kartaltepe et al 2015)` == "yes" & `Mostly Spheroid` == "no" & `Mostly Irregular` == "no"')
    disk_irr = tab.query('`Mostly Disk (4.3 Kartaltepe et al 2015)` == "yes" & `Mostly Spheroid` == "no" & `Mostly Irregular` == "yes"')
    disk_sph = tab.query('`Mostly Disk (4.3 Kartaltepe et al 2015)` == "yes" & `Mostly Spheroid` == "yes" & `Mostly Irregular` == "no"')
    spheroid_only = tab.query('`Mostly Disk (4.3 Kartaltepe et al 2015)` == "no" & `Mostly Spheroid` == "yes" & `Mostly Irregular` == "no"')
    irr_sph = tab.query('`Mostly Disk (4.3 Kartaltepe et al 2015)` == "no" & `Mostly Spheroid` == "yes" & `Mostly Irregular` == "yes"')
    irregular_only = tab.query('`Mostly Disk (4.3 Kartaltepe et al 2015)` == "no" & `Mostly Spheroid` == "no" & `Mostly Irregular` == "yes"')
    disk_irr_sph = tab.query('`Mostly Disk (4.3 Kartaltepe et al 2015)` == "yes" & `Mostly Spheroid` == "yes" & `Mostly Irregular` == "yes"')
    orphans = tab.query('`Mostly Disk (4.3 Kartaltepe et al 2015)` == "no" & `Mostly Spheroid` == "no" & `Mostly Irregular` == "no"')

    ### add them to general array
    morph_types = np.array([disk_only,spheroid_only,irregular_only,disk_sph,disk_irr,irr_sph,disk_irr_sph,orphans],dtype=object)
    labels = ['Mostly Disk', 'Mostly Spheroid', 'Mostly Irregular','Disk+Spheroid','Disk+Irregular','Irregular+Spheroid','Disk+Irregular+Spheroid','Unclassifiable']
    
    classes = np.full(table['Properties'].shape[0],fill_value='None',dtype=object)
    z = 0
    for i in morph_types:
        ind = i.index.tolist()
        classes[ind] = labels[z]
        z+=1

    tab['Classification'] = classes

    return morph_types,labels
morph_types, labels = add_classifications(table)

##### get the merger variants for each classification
##### for now we're including the 'maybe' flag, but in the future we should exclude it (or just exclude if from the sheet)
v = 0
mergers = np.empty(len(morph_types),dtype=object)
for i in morph_types:
    mergers[v] = i.query('`Merger (flag threshold)` == "yes" | `Merger (flag threshold)` == "maybe"')
    v+=1


#### plotting arrays
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
print(f'Total orphan sources: {orphans.shape[0]} ({orphans["Properties"]})')

#### show # in each category
c=0
for i in morph_types:
    print(f'{labels[c]}: {i.shape[0]}')
    c+=1

morph_agn = [i[['AGN(%)']] for i in morph_types]
morph_z = [i[['z']] for i in morph_types]  
morph_clumps = [i.query('`Has Clumps (flag)`=="yes"') for i in morph_types]  
morph_spiral = [i.query('`Spiral Arms (flag)`=="yes"') for i in morph_types]

def agn_frac_hist():
    ##### Make a plot of galaxy type by agn fraction
    histbins = 5
    plt.title("AGN Fraction by Galaxy Type (Disk-only w/ Clumps)")
    plt.ylabel('Count')
    plt.xlabel('AGN Fraction (%)')
    bins=np.histogram(np.hstack((disk_agn,sph_agn,irr_agn,disk_sph_agn,disk_irr_agn,disk_irr_sph_agn,orphans['AGN(%)'])), bins=histbins)[1] #get the bin edges
    # plt.hist(disk_agn, label='Mostly Disk',bins=bins,ec='blue',lw=1,fill=False,histtype='step')
    # plt.hist(sph_agn, label='Mostly Spheroid',bins=bins,ec='red',lw=1,fill=False)
    # plt.hist(irr_agn, label='Mostly Irregular',bins=bins,ec='green',lw=1,fill=False)
    # plt.hist(disk_sph_agn, label='Disk+Spheroid',bins=bins,ec='magenta',lw=1,fill=False)
    # plt.hist(disk_irr_agn, label='Disk+Irregular',bins=bins,ec='orange',lw=1,fill=False)
    # plt.hist(irr_sph_agn, label='Irregular+Spheroid',bins=bins,ec='yellow',lw=1,fill=False)
    # plt.hist(disk_irr_sph_agn, label='Disk+Spheroid+Irregular',bins=bins,ec='cyan',lw=1,fill=False)
    # plt.hist(orphans['AGN(%)'],label='Unclassifiable',bins=bins,ec='gray',lw=1,fill=False)
    # plt.legend()
    # plt.vlines([20,80])
    # agn_filename = 'agn_frac_spheroid_only'
    # plt.savefig(f'{working_directory}{agn_filename}',dpi = 300)
    # plt.show()

    ## plotting all of the morph_clumps in array
    # plot the fraction of disks that have clumps to histogram
    c=0
    for i in morph_clumps:
        plt.ylabel('Count')
        plt.xlabel('AGN Fraction (%)')
        plt.hist(morph_agn[c], label=labels[c],bins=bins,lw=1,fill=False,histtype='step')
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
        plt.hist(morph_agn[c], label=labels[c],bins=bins,lw=1,fill=False,histtype='step')
        plt.hist(i['AGN(%)'],bins=bins,fill=True,label='Has Arms')
        plt.title(f"AGN Fraction by Galaxy Type ({labels[c]} w/ Spiral Arms)")
        plt.legend()
        plt.savefig(f'{working_directory}agn-hist-w-spirals-{labels[c]}',dpi = 300)
        # plt.show()
        plt.close()
        c+=1


    # ### get the names and display cutouts of the jades rgb images for different ranges of each
    # irr = irregular_only.sort_values(by='AGN(%)')
    # print('IRREGULARS -----')
    # print(irr[['Properties','AGN(%)']]) # get both properties & agn frac
    # sph = spheroid_only.sort_values(by='AGN(%)')
    # print('SPHEROIDS -----')
    # print(sph[['Properties','AGN(%)']]) # get both properties & agn frac
    # disk = disk_only.sort_values(by='AGN(%)')
    # print('DISKS -----')
    # print(disk[['Properties','AGN(%)']]) # get both properties & agn frac


def type_z_hist():
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


def agn_frac_z():
    #### make a scatter plot of the agn fraction vs redshift to showcase if there is a relation there?
    
    
    ### Calculate the % morphology in each agn section (<20%, 20<x<80%, >80%)
    
    ### old but save for later bc this shorthand is super OP
    # sf_sample = [i.query('`AGN(%)`<=20') for i in morph_types]
    # comp_sample = [i.query('`AGN(%)`>20 & `AGN(%)`<80') for i in morph_types]
    # pure_agn_sample = [i.query('`AGN(%)`>=80') for i in morph_types]
    bin_colors = {
        'disk only': '#1f77b4',
        'spheroid only': '#ff7f0e',
        'irregular only': '#2ca02c',
        'disk+spheroid': '#d62728',
        'disk+irregular': '#9467bd',
        'spheroid+irregular': '#8c564b',
        'disk+spheroid+irregular': '#e377c2',
        'Unclassifiable': '#7f7f7f'
    }
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
    sf_sample_table = table.query('`AGN(%)`<=20')
    comp_sample_table = table.query('`AGN(%)`>20 & `AGN(%)`<80')
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
    hist_bins = np.histogram(np.hstack((morph_types[0]['AGN(%)'],morph_types[1]['AGN(%)'],morph_types[2]['AGN(%)'],morph_types[3]['AGN(%)'])),bins=histbins)[1] # ideally would put all them but I'm lazy and it looks good as is
    for i in mergers:
        plt.hist(i['AGN(%)'],bins=hist_bins,lw=1.5,fill=False,histtype='step')
    plt.legend(labels,loc='best',prop={'size': 8})
    plt.xlabel('AGN(%)')
    plt.ylabel('Count')
    plt.title('AGN(%) of Mergers')
    plt.savefig(f'{working_directory}agn_merger_hist',dpi=300)
    plt.show()

def agn_frac_merger_scatter():
    #### make a scatter plot of the agn fraction vs merger (% agreement) with color coded classifications
    for i in mergers:
        plt.scatter(i['AGN(%)'],i['Merger? (%)'])
    plt.xlabel('AGN(%)')
    plt.ylabel('Agrement on Merger (%)')
    plt.xlim([0,50])
    plt.ylim([0,100])
    plt.axhline(66,linestyle='--',color='darkgray')
    plt.axhline(66/2,linestyle='--',color='lightgray')
    labels.append('Yes Threshold')
    labels.append('Maybe Threshold')
    plt.title('AGN(%) of Mergers')
    plt.legend(labels,loc='upper right',prop={'size': 7})
    plt.savefig(f'{working_directory}agn_merger_scatter',dpi=300)
    plt.show()
    labels.remove('Yes Threshold')
    labels.remove('Maybe Threshold')

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

    '''
    sub_sil = df.query('`Sub-sample`=="Sil AGN"')[['ID','AGN^a']]
    sub_feat = df.query('`Sub-sample`=="Feat AGN"')[['ID','AGN^a']]
    sub_unknown = df.query('`Sub-sample`=="..."')[['ID','AGN^a']]

    # frame = table.where(table['Properties']==sub_sil)
    # frame = table[(table.Properties == sub_sil)]
    sub_sil_classify = table[(table['Properties'].isin(sub_sil['ID']))].sort_values(by='Properties',ascending=True)
    sub_feat_classify = table[(table['Properties'].isin(sub_feat['ID']))].sort_values(by='Properties',ascending=True)
    sub_unknown_classify = table[(table['Properties'].isin(sub_feat['ID']))].sort_values(by='Properties',ascending=True)

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

nir_filter = 'f444w'
stat_measures = f"research/statmorph_output/{nir_filter}_6_JADES_statmorph_measurements.tsv"
stats = pd.read_csv(stat_measures,sep='\t')
stats.sort_values(by='ID',ascending=True,inplace=True)
agn_dominated = table.query('`AGN(%)`>=50')
sf_dominated = table.query('`AGN(%)`<50')

def A_S_scatter():
    ### scatter the asymmetry vs smoothness (clumpiness) measurements
    # agn_stats= stats[(stats['ID'].str.contains(id,regex=False))].iloc[0]
    agn_stats = stats[(stats['ID'].isin(agn_dominated['Properties']))]
    sf_stats = stats[stats['ID'].isin(sf_dominated['Properties'])]

    plt.scatter(agn_stats['Smoothness (S)'],agn_stats['Asymmetry (A)'],label='AGN')
    plt.scatter(sf_stats['Smoothness (S)'],sf_stats['Asymmetry (A)'],label='SF')
    plt.xlabel('Clumpiness (S)')
    plt.ylabel('Asymmetry (A)')
    plt.legend()
    plt.ylim((-.1,.9))
    plt.xlim((-.025,.125))
    plt.title(f'{nir_filter.upper()} A vs S')
    plt.gca().invert_yaxis()
    plt.show()



# agn_frac_hist()
# type_z_hist()

# agn_frac_z()
# testing_agn_frac_z()

# agn_frac_merger_hist()
# agn_frac_merger_scatter()
# compare_subsample()
A_S_scatter()




disclaimer()