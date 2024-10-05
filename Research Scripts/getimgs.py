#### imports for making galaxy img cutouts
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd

#### imports for making galaxy img cutouts
import astropy.units as units
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.nddata import Cutout2D
from astropy.io import fits
from photutils.segmentation import SegmentationImage

from import_data import get_GS,get_GN,get_sources

### getting images for the poster/large cutout image on paper
def get_grizli_images(filter=None):
    ######## store everything in respective array location for the following Nircam filters
    if(filter!=None):
        nir_filters = [filter]
    else:
        nir_filters = ['f150w','f200w','f277w','f356w','f444w']


    gs_scis = np.zeros(len(nir_filters),dtype=object)
    gs_wcs_coords = np.zeros(len(nir_filters),dtype=object)
    gs_exp_wcs_coords = np.zeros(len(nir_filters),dtype=object)
    gs_wht_wcs = np.zeros(len(nir_filters),dtype=object)
    gs_whtmaps = np.zeros(len(nir_filters),dtype=object)
    gs_expmaps = np.zeros(len(nir_filters),dtype=object)

    gs_noise_vars = np.zeros((len(nir_filters),4))  ### layout is PHOTMJSR, PHOTSCAL, PHOTFNU, OPHOTFNU
    gn_noise_vars = np.zeros((len(nir_filters),4))  ### layout is PHOTMJSR, PHOTSCAL, PHOTFNU, OPHOTFNU
    

    gn_scis = np.zeros(len(nir_filters),dtype=object)
    gn_wcs_coords = np.zeros(len(nir_filters),dtype=object)
    gn_exp_wcs_coords = np.zeros(len(nir_filters),dtype=object)
    gn_wht_wcs = np.zeros(len(nir_filters),dtype=object)
    gn_whtmaps = np.zeros(len(nir_filters),dtype=object)
    gn_expmaps = np.zeros(len(nir_filters),dtype=object)
    
    gs_jades_seg = None
    full_jades_seg = None # for testing with grizli so that we use thE JADES SEGMAP
    gs_jades_wcs = None #   ^^^^^^^^^^^
    ### try to use the grizli data reduction otherwise use jades

    ind=0
    for f in nir_filters:
        hdul = fits.open(f'/home/robbler/research/grizli-reductions/gds-grizli-v7.2-{f}-clear_drc_sci.fits')
        gs_whtmap_general = fits.open(f'/home/robbler/research/grizli-reductions/gds-grizli-v7.2-{f}-clear_drc_wht.fits')
        gs_whtmaps[ind] = gs_whtmap_general[0].data
        gs_wht_wcs[ind] = WCS(gs_whtmap_general[0].header)
        gs_expmap_general = fits.open(f'/home/robbler/research/grizli-reductions/gds-grizli-v7.2-{f}-clear_drc_exp.fits')
        gs_expmaps[ind] = gs_expmap_general[0].data
        gs_exp_wcs_coords[ind] = WCS(gs_expmap_general[0].header)
        gs_expmap_general.close()
        gs_whtmap_general.close()
        # gs_jades_seg = fits.getdata('/home/robbler/research/grizli-reductions/gds-grizli-v7.2-ir_seg.fits')
        gs_jades_seg = fits.getdata('/home/robbler/research/JADES_Fits_Maps/GS/hlsp_jades_jwst_nircam_goods-s-deep_segmentation_v2.0_drz.fits')
        full_jades_seg = fits.open('/home/robbler/research/JADES_Fits_Maps/GS/hlsp_jades_jwst_nircam_goods-s-deep_segmentation_v2.0_drz.fits')
        
        gs_scis[ind] = hdul[0].data


        #### adding poisson noise parameters
        gs_header = hdul[0].header

        gs_noise_vars[ind][0] = gs_header['PHOTMJSR']
        gs_noise_vars[ind][1] = gs_header['PHOTSCAL']
        gs_noise_vars[ind][2] = gs_header['PHOTFNU']
        gs_noise_vars[ind][3] = gs_header['OPHOTFNU']

        gs_wcs_coords[ind] = WCS(hdul[0].header)
        ###### for testing with grizli sources so that we use the JADES SEGMAP instead
        gs_jades_wcs = WCS(full_jades_seg[0].header)
        
        # gn_wcs_coords = WCS(hdul[0].header)
        full_jades_seg.close()

        # catalog file
        tmp = fits.open('/home/robbler/research/JADES catalog/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits')

        # generate PSF for current filter
        # nc = webbpsf.NIRCam()
        # nc.filter = 'F444W'
        # psf = nc.calc_psf(oversample=4,display=True)
        # print(psf)
        # exit()

        hdul.close()

        # catalog file
        jades_gs_sizes = Table(tmp['SIZE'].data)
        gs_jades_catalog = Table(tmp['FLAG'].data) # extension #2 "FLAG" data
        tmp.close()

        # print(gn_whtmap[0].header['FLAGS_WEIGHT'])
        # exit()

        #################### open goods north files ... ###################################################################################
        # starting with f444 because it should have the clearest source (maybe not clear detail though)
        #hdul = fits.open('JADES_Fits_Maps/hlsp_jades_jwst_nircam_goods-s-deep_f444w_v2.0_drz.fits')
        
        ### try to use the grizli data reduction otherwise use jades
        gn_jades_seg = None
        gn_jades_wcs = None

        hdul = fits.open(f'/home/robbler/research/grizli-reductions/gdn-grizli-v7.3-{f}-clear_drc_sci.fits')
        gn_whtmap_general = fits.open(f'/home/robbler/research/grizli-reductions/gdn-grizli-v7.3-{f}-clear_drc_wht.fits')
        gn_whtmaps[ind] = gn_whtmap_general[0].data
        gn_wht_wcs[ind] = WCS(gn_whtmap_general[0].header)
        gn_expmap_general = fits.open(f'/home/robbler/research/grizli-reductions/gdn-grizli-v7.3-{f}-clear_drc_wht.fits')
        gn_expmaps[ind] = gn_expmap_general[0].data
        gn_exp_wcs_coords[ind] = WCS(gn_expmap_general[0].header)

        gn_whtmap_general.close()
        gn_expmap_general.close()
        gn_jades_seg = fits.getdata('research/JADES_Fits_Maps/GN/hlsp_jades_jwst_nircam_goods-n_segmentation_v1.0_drz.fits')
        full_jades_seg = fits.open('research/JADES_Fits_Maps/GN/hlsp_jades_jwst_nircam_goods-n_segmentation_v1.0_drz.fits')

        gn_scis[ind] = hdul[0].data


        ### adding poisson noise variables
        gn_header = hdul[0].header

        gn_noise_vars[ind][0] = gn_header['PHOTMJSR']
        gn_noise_vars[ind][1] = gn_header['PHOTSCAL']
        gn_noise_vars[ind][2] = gn_header['PHOTFNU']
        gn_noise_vars[ind][3] = gn_header['OPHOTFNU']

        gn_wcs_coords[ind] = WCS(hdul[0].header)
        ###### for testing with grizli sources so that we use the JADES SEGMAP instead
        gn_jades_wcs = WCS(full_jades_seg[0].header)
        # gn_wcs_coords = WCS(hdul[0].header)
        full_jades_seg.close()

        # catalog file
        tmp = fits.open('research/JADES catalog/hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits')

        # generate PSF for current filter
        # nc = webbpsf.NIRCam()
        # nc.filter = 'F444W'
        # psf = nc.calc_psf(oversample=4,display=True)
        # print(psf)
        # exit()

        hdul.close()
        # catalog file
        jades_gn_sizes = Table(tmp['SIZE'].data)
        gn_jades_catalog = Table(tmp['FLAG'].data) # extension #2 "FLAG" data
        tmp.close()

        ind+=1
    sources = get_sources(full=True)


    #### loop through the filters here or something
    for i in nir_filters:
        f_s = sources.query(f'`Obs Filter`=="{i}"')
        filter_sources = f_s.sort_values(by='AGN^a')
        # f277_sources = sources.query('`Obs Filter`=="f277w"')
        # f277_sources.sort_values(by='AGN^a',inplace=True)

        search_ids = filter_sources['ID']


        images = np.zeros(len(search_ids),dtype=object)
        whts = np.zeros(len(search_ids),dtype=object)
        segs = np.zeros(len(search_ids),dtype=object)

        ids = np.zeros(len(search_ids),dtype=str)
        agns = np.zeros(len(search_ids),dtype=float)

        ind=axind=0
        # try to make a square grid based on the length of the array
        square = round(np.ceil(np.sqrt(len(filter_sources))))
        nrows, ncols = square,square
        # new_shape = (nrows,ncols)
        fig, axs = plt.subplots(nrows,ncols, figsize=(10, 10))
        axs=axs.flatten()
        for id in search_ids:
            axind=ind
            row = []

            data_source = sources[(sources['ID'].str.contains(id,regex=False))].iloc[0]
            source_row = sources.loc[sources['ID'] == id]
            source_filter = source_row['Obs Filter'].iloc[0]
            index = nir_filters.index(source_filter)

            size = 5 ### constant size, this should be in arc seconds, but unsure if the conversion is correct



            search_id = id
            jades_catalog,sci,jades_seg,wcs_coords,jades_wcs_coords,whtmap = (None for i in range(6))
            if(id.startswith('GN')):
                jades_catalog = gn_jades_catalog
                sci = gn_scis[index]
                wcs_coords = gn_wcs_coords[index]
                whtmap = gn_whtmaps[index]
                wht_wcs = gn_wht_wcs[index]
                jades_seg = gn_jades_seg
                jades_wcs_coords = gn_jades_wcs


            if(id.startswith('GS')):
                jades_catalog = gs_jades_catalog
                sci = gs_scis[index]
                wcs_coords = gs_wcs_coords[index]
                whtmap = gs_whtmaps[index]
                wht_wcs = gs_wht_wcs[index]
                jades_seg = gs_jades_seg
                jades_wcs_coords = gs_jades_wcs



            jades_id_raw = data_source['JADES ID']
            jades_id = 0
            if(' ' in jades_id_raw):
                jades_id = int(jades_id_raw.split(' ')[0]) # splitting because sometimes I have a secondary source on there, but the main is always listed first (this only happens 2 times in our sample)
            else:
                jades_id = int(jades_id_raw)

            source = jades_catalog[jades_catalog['ID']==jades_id]
            position = SkyCoord(source['RA'],source['DEC'],unit='deg')

            if(id=='GS_IRS37'):
                ## skip this source because it's not in the map
                continue
            images[ind] = Cutout2D(sci,position,size*units.arcsec,wcs=wcs_coords,copy=True).data
            whts[ind] = Cutout2D(whtmap,position,size*units.arcsec,wcs=wht_wcs,copy=True).data
            if(i=='f150w' or i=='f200w'):
                # different scaling for segmentation map on 150 & 200
                segs[ind] = Cutout2D(jades_seg,position,1.5*size*units.arcsec,wcs=jades_wcs_coords,copy=True).data
            else:
                segs[ind] = Cutout2D(jades_seg,position,0.75*size*units.arcsec,wcs=jades_wcs_coords,copy=True).data

            # seg_img = SegmentationImage(segs[ind])
            # segs[ind] = seg_img
            # plt.imshow(segs[ind],origin='lower')
            # plt.show()
            # exit()

            # segs[ind] = Cutout2D(jades_seg,position,size*units.arcsec,wcs=jades_wcs_coords,copy=True).data
            ids[ind] = id
            agns[ind] = source_row['AGN^a']

            
            # rms = np.sqrt(np.mean(np.sqrt(1/whts[ind])))
            # c = np.divide(1, whts[ind], out=np.zeros_like(a), where=b!=0)
            # noise_arr = np.where(np.sqrt(1/whts[ind])!=np.inf)
            noise_arr = np.sqrt(1/whts[ind])
            noise_arr[noise_arr==np.inf] = 0
            rms = np.sqrt(np.mean(noise_arr))
            # print(rms)



            if(ind>19 and i=='f277w'):
                ## rigging the plot to show the last 2 agn sources of f277 on the right instead of left
                axind+=2

            # print(id)
            # plt.imshow(segs[ind],origin='lower')
            # plt.show()
            # exit()
            axs[axind].imshow(images[ind],origin='lower',norm=colors.PowerNorm(gamma=0.95,vmin=-3*rms,vmax=35*rms),cmap='gray_r')
            # axs[axind].imshow(segs[ind],origin='lower')
            # axs[axind].contour(segs[ind], levels=[jades_id-5,jades_id,jades_id+5], colors='red', linewidths=1)
            axs[axind].contour(segs[ind], levels=[jades_id-1,jades_id,jades_id+1], colors=['red'], linewidths=1)
            # axs[axind].contour(segs[ind], levels=[jades_id], colors=['red'], alpha=1,linewidths=1)
            # axs[axind].contour(segs[ind], levels=[jades_id-50,jades_id+50],colors=['blue'], alpha=1, linewidths=1)


            # might need to adjust the scaling for each of them, also need to make sure the size makes sense for each too
            axs[axind].text(5,5, id, color='black', fontsize=12, fontweight='bold')

            # axs[ind].text(5,55,agns[ind],color='black',fontsize=12,fontweight='bold')
            for s in axs[axind].spines.values():
                if(agns[ind]>=80):
                    s.set_edgecolor('blue')
                elif(agns[ind]>20 and agns[ind]<80):
                    s.set_edgecolor('orange')    
                elif(agns[ind]<=20):
                    s.set_edgecolor('green')
                s.set_linewidth(2)

            axs[axind].set_xticks([])
            axs[axind].set_yticks([])

            
            # spine = axs[ind].spines.values()
            # spine.set_edgecolor('blue')
            # axs[ind].axis('off')
            # plt.tight_layout(pad=1)

            ##### maybe something off with the scaling of the gs_irs12 image?
            ind+=1

        for a in range(len(axs)):
            if not axs[a].has_data():
                fig.delaxes(axs[a])
            
        # for i in range(len(search_ids),len(axs)):
        #     fig.delaxes(axs[i])

        fig.suptitle(f'{i} Cutouts by AGN Fraction ({size}")',size=16)
        plt.tight_layout(pad=1)
        plt.show()
        ### save the plt here
        plt.close()

    return
        
get_grizli_images(filter='f277w')



################# below this is depreciated #################################


### getting images for the visual form
def get_GS_images(filter='f150w',name=None):
    '''
    Function to set the current field to GS before running the cutout loop through each source.
    Optional parameter: Filter -> like 'f150w', 'f277w', etc.
    Optional Parameter: Name -> source name to search for within the sample (ex. GS_IRS81)
    Returns the filter name (as inputted)
    '''
    filename = (f'/home/robbler/research/JADES_Fits_Maps/GS/hlsp_jades_jwst_nircam_goods-s-deep_{filter}_v2.0_drz.fits')
    if(os.path.exists(filename)):
        with fits.open(filename) as hdul:
        #hdul = fits.open('/home/robbler/research/JADES_Fits_Maps/GS/hlsp_jades_jwst_nircam_goods-s-deep_f150w_v2.0_drz.fits')

        # .fits file (image version)
            fits_header = hdul[1].header # need this without the wcs for the pyregion library
            sci = hdul[1].data # for getting the science image
            wcs_coords = WCS(hdul[1].header) # getting wcs from header
        #hdul.close()
    else:

        raise Exception(f"Filter provided was invalid: '{filter}'")


    # catalog file
    tmp = fits.open('/home/robbler/research/JADES catalog/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits')
    jades_catalog = Table(tmp['FLAG'].data) # extension #2 "FLAG" data
    tmp.close()

    source_names, jades_names, source_coords = get_GS()

    regions_file = '/home/robbler/research/Galaxies/Locations/GS_regions_file.reg'
    ### Specific filter settings
    ## f150w
    #im = ax.imshow(img,origin='lower',vmin=(m-s),vmax=(m+5*s),cmap='gray_r')
    ## f277w
    #im = ax.imshow(img,origin='lower',norm=colors.PowerNorm(gamma=0.5,vmin=(m-s),vmax=(m+5*s)),cmap='gray_r')

    ###################
    """
    Now take the cutout images for all sources (from specified sources array)
    """
    ###################

    c=0 # for testing normalization

    ### only show the cutout for the named galaxy if name exists
    if(name!=None):
        source_names = np.array(source_names)
        ind = np.where(source_names==str(name))[0][0]
        # cut down the names and coords to just the desired source coords
        source_names=np.array([name])
        source_coords = np.array([[source_coords[ind,0],source_coords[ind,1]]])
        #source_coords = np.array([[53.15216667,-27.77519444]])
        
    
    # detection radius ~ 1" for Kirkpatrick data, so we want anything within a ~3" box

    for i in source_names:
        size = 3

        #source = jades_catalog[jades_catalog['ID']==id]
        #position = SkyCoord(source['RA'],source['DEC'],unit='deg')
        position = SkyCoord(source_coords[c,0],source_coords[c,1],unit='deg')
        # position = SkyCoord(53.15219177,-27.77526654,unit='deg')
        image_base = Cutout2D(sci,position,size*units.arcsec,wcs=wcs_coords)
        img = image_base.data

        #fig, ax = plt.subplots(figsize=(6,6),subplot_kw={'projection':wcs_coords})
        fig = plt.figure(figsize=(8,8))
        #print(image_base.wcs)
        ax = fig.add_subplot(1,1,1, projection=image_base.wcs)
        #fig = plt.figure()
        #ax = WCSAxes(fig,[0.1, 0.1, 0.8, 0.8], wcs=wcs_coords)
        #fig.add_axes(ax)


        ## practice making edits to fits files in astropy
        ## https://learn.astropy.org/tutorials/FITS-images.html

        #from matplotlib.colors import LogNorm
        # display the image to see what we're working with
        m = np.mean(img)
        s = np.std(img)


        if(filter=='f150w'):
            #im = ax.imshow(img,origin='lower',vmin=(m-s),vmax=(m+5*s),cmap='gray_r')
            #im = ax.imshow(img,origin='lower',vmin=(m-s),vmax=(m+8*s),cmap='gray_r') (m-0.5*s)
            im = ax.imshow(img,origin='lower',norm=colors.PowerNorm(gamma=0.5,vmin=(m-s),vmax=(m+8*s)),cmap='gray_r')
        elif(filter=='f277w'):
            # just so happens to be the same filter settings are the best?
            im = ax.imshow(img,origin='lower',norm=colors.PowerNorm(gamma=0.5,vmin=(m-s),vmax=(m+8*s)),cmap='gray_r')
        else:
            im = ax.imshow(img,origin='lower',vmin=(m-s),vmax=(m+5*s),cmap='gray_r')
            raise Exception(f'No valid filter parameters found, please adjust filter. {filter=}')
        # this is for the 150 filter
        #im = ax.imshow(img,origin='lower',vmin=(m-s),vmax=(m+5*s),cmap='gray_r')
        # this might be for 277?
        #im = ax.imshow(img,origin='lower',norm=colors.PowerNorm(gamma=0.5,vmin=(m-s),vmax=(m+5*s)),cmap='gray_r')

        plt.title(f'{source_names[c]} ({filter})') # this should be i inside the loop, assuming the array you're looping through follows this structure
    
        # a little computationally expensive to loop through all regions for each source, but it seems like too much work to fix it
        ##### label the filter with the regions file made from "read_jades_catalog.py"
        regions = Regions.read(regions_file, format='ds9')

        for i, reg in enumerate(regions):
            pix_reg = reg.to_pixel(image_base.wcs)       # For plotting, Region objects have to be in pixel coordinates
            pix_reg.visual['linewidth'] = 1
            pix_reg.plot(ax=ax)
        plt.colorbar(im,ax=ax)
        figname = (f'{source_names[c]}_{filter}')

        if(name!=None):
            plt.show()
            # plt.savefig(f'/home/robbler/research/Galaxies/cutouts/GS/testing_{figname}')
        else:
            plt.savefig(f'/home/robbler/research/Galaxies/cutouts/GS/{figname}')
            print(f'Saved figure: {figname}')
        c+=1

# get_GS_images('f150w',name='GS_IRS81')

### getting images for the visual form
def get_GN_images(filter='f150w',name=None):
    '''
    Function to set the current field to GN before running the cutout loop through each source.
    Optional parameter: Filter -> like 'f150w', 'f277w', etc.
    Optional parameter: Name -> ID of source you want to look at like 'GN_IRS11'
    Returns the filter name (as inputted)
    '''
    filename = (f'/home/robbler/research/JADES_Fits_Maps/GN/hlsp_jades_jwst_nircam_goods-n_{filter}_v1.0_drz.fits')
    if(os.path.exists(filename)):
        with fits.open(filename) as hdul:
        #hdul = fits.open('/home/robbler/research/JADES_Fits_Maps/GS/hlsp_jades_jwst_nircam_goods-s-deep_f150w_v2.0_drz.fits')

        # .fits file (image version)
            fits_header = hdul[1].header # need this without the wcs for the pyregion library
            sci = hdul[1].data # for getting the science image
            wcs_coords = WCS(hdul[1].header) # getting wcs from header
        #hdul.close()
    else:

        raise Exception(f"Filter provided was invalid: '{filter}'")


    # catalog file
    tmp = fits.open('/home/robbler/research/JADES catalog/hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits')
    jades_catalog = Table(tmp['FLAG'].data) # extension #2 "FLAG" data
    tmp.close()
     
    source_names, jades_name_data, source_coords = get_GN()
    # this formatting is different than GS [:,0] here vs [0,:] before (for finding RA as ex.)
    #source_coords = get_ints(coord_data)
    #source_names = get_strings(name_data)
    regions_file = '/home/robbler/research/Galaxies/Locations/GN_regions_file.reg'
    ### Specific filter settings (notes for later?)
    ## f150w
    #im = ax.imshow(img,origin='lower',vmin=(m-s),vmax=(m+5*s),cmap='gray_r')
    ## f277w
    #im = ax.imshow(img,origin='lower',norm=colors.PowerNorm(gamma=0.5,vmin=(m-s),vmax=(m+5*s)),cmap='gray_r')

    ###################
    """
    Now take the cutout images for all sources (from specified sources array)
    """
    ###################
    
    
    c=0 # for testing normalization

    ### if name is set then only get and display the image of this single galaxy, otherwise do all them
    # name = 'GN_IRS31^e'
    if(name!=None):
        source_names = np.array(source_names)
        ind = np.where(source_names==str(name))[0][0]
        # cut down the names and coords to just the desired source coords
        source_names=np.array([name])
        source_coords = np.array([[source_coords[ind,0],source_coords[ind,1]]])


    # detection radius ~ 1" for Kirkpatrick data, so we want anything within a ~3" box
    for i in source_names:
        size = 3

        #source = jades_catalog[jades_catalog['ID']==id]
        #position = SkyCoord(source['RA'],source['DEC'],unit='deg')
        position = SkyCoord(source_coords[c,0],source_coords[c,1],unit='deg')
        image_base = Cutout2D(sci,position,size*units.arcsec,wcs=wcs_coords)
        img = image_base.data

        #fig, ax = plt.subplots(figsize=(6,6),subplot_kw={'projection':wcs_coords})
        fig = plt.figure(figsize=(8,8))
        #print(image_base.wcs)
        ax = fig.add_subplot(1,1,1, projection=image_base.wcs)
        #fig = plt.figure()
        #ax = WCSAxes(fig,[0.1, 0.1, 0.8, 0.8], wcs=wcs_coords)
        #fig.add_axes(ax)


        ## practice making edits to fits files in astropy
        ## https://learn.astropy.org/tutorials/FITS-images.html

        #from matplotlib.colors import LogNorm
        # display the image to see what we're working with
        m = np.mean(img)
        s = np.std(img)


        if(filter=='f150w'):
            #im = ax.imshow(img,origin='lower',vmin=(m-s),vmax=(m+5*s),cmap='gray_r')
            #im = ax.imshow(img,origin='lower',vmin=(m-s),vmax=(m+8*s),cmap='gray_r') (m-0.5*s)
            im = ax.imshow(img,origin='lower',norm=colors.PowerNorm(gamma=0.5,vmin=(m-s),vmax=(m+8*s)),cmap='gray_r')
            # im = ax.imshow(img,origin='lower',norm=colors.PowerNorm(gamma=1.,vmin=(m-s),vmax=(m+2*s)),cmap='gray_r')
        elif(filter=='f277w'):
            # just so happens to be the same filter settings are the best?
            im = ax.imshow(img,origin='lower',norm=colors.PowerNorm(gamma=0.5,vmin=(m-s),vmax=(m+8*s)),cmap='gray_r')
            # im = ax.imshow(img,origin='lower',norm=colors.PowerNorm(gamma=1.,vmin=(m-s),vmax=(m+2*s)),cmap='gray_r')
        else:
            im = ax.imshow(img,origin='lower',vmin=(m-s),vmax=(m+5*s),cmap='gray_r')
            raise Exception(f'No valid filter parameters found, please adjust filter. {filter=}')
        # this is for the 150 filter
        #im = ax.imshow(img,origin='lower',vmin=(m-s),vmax=(m+5*s),cmap='gray_r')
        # this might be for 277?
        #im = ax.imshow(img,origin='lower',norm=colors.PowerNorm(gamma=0.5,vmin=(m-s),vmax=(m+5*s)),cmap='gray_r')

        plt.title(f'{source_names[c]} ({filter})') # this should be i inside the loop, assuming the array you're looping through follows this structure
        ## [[kirk id 1, kirk id 2, etc],[jades id 1, jades id 2, etc]]
        ## it would work this easily because we have all the kirk and jades sources matched to only one source for all but a few so that should work fine
    
        # a little computationally expensive to loop through all regions for each source, but it seems like too much work to fix it
        ##### label the filter with the regions file made from "read_jades_catalog.py"
        regions = Regions.read(regions_file, format='ds9')

        for i, reg in enumerate(regions):
            pix_reg = reg.to_pixel(image_base.wcs)       # For plotting, Region objects have to be in pixel coordinates
            pix_reg.visual['linewidth'] = 1
            pix_reg.plot(ax=ax)
        plt.colorbar(im,ax=ax)
        figname = (f'{source_names[c]}_{filter}')
        if(name!=None):
            plt.show()
            # plt.savefig(f'/home/robbler/research/Galaxies/cutouts/GN/testing_{figname}')
        else:
            plt.savefig(f'/home/robbler/research/Galaxies/cutouts/GN/{figname}')
        print(f'Saved figure: {figname}')
        plt.close()
        c+=1


# get_GN_images(filter='f277w')