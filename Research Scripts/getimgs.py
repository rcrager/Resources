#### imports for making galaxy img cutouts
import astropy.units as units
from astropy.wcs import WCS
from astropy.coordinates import Angle, SkyCoord
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from regions import Regions
import os
from import_data import get_GN, get_GS

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

get_GS_images('f150w',name='GS_IRS81')


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