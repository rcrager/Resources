import statmorph

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.visualization import simple_norm
from astropy.modeling.models import Sersic2D
from astropy.convolution import convolve, Gaussian2DKernel
import pandas as pd

#### imports for making galaxy img cutouts
import astropy.units as units
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, join
from astropy.nddata import Cutout2D
from astropy.io import fits,ascii
from astropy.stats import sigma_clipped_stats
import warnings

from import_data import get_GS,get_GN,get_sources,disclaimer

### imports for making the segmentation map (if need be)
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources
from photutils.background import Background2D, MedianBackground
from statmorph.utils.image_diagnostics import make_figure
from photutils.segmentation import make_2dgaussian_kernel, SegmentationImage
from astropy.stats import SigmaClip
#sigma_clipping = SigmaClip(sigma=sigma)

# for finding the psf
import webbpsf

from astropy.utils.exceptions import AstropyWarning
#warnings.simplefilter('ignore', category=AstropyWarning)
# following https://statmorph.readthedocs.io/en/latest/notebooks/tutorial.html
## to learn it

# import sys                          # <-comment out to see the output
# tmp = sys.stdout                    # <-comment out to see the output

# output_file = '/home/robbler/research/statmorph_output/log_statmorph_output.txt'
# sys.stdout = open(output_file, "w")   # <-comment out to see the output


def run_in_jades(filter='f277w',output_name='JADES_statmorph_measurements',use_grizli=False):
    ###############################
    """
    Get statmorph fit for a jades source
    Input:  filter: NIRCam Filter to select
            output_name: Filename to output the table of measurements to (saved filename is inclusive of filter name)
    """
    ###############################
    ## making a sample image that is going to be analyzed (with sersic 2d)
    '''
    ny, nx = 240,240
    y,x = np.mgrid[0:ny, 0:nx]
    sersic = Sersic2D(amplitude=1, r_eff=20, n=2.5, x_0=120.5, y_0=96.5, ellip=0.5, theta=0.5)

    image = sersic(x,y)
    #plt.imshow(image, cmap='gray', origin='lower', norm=simple_norm(image,stretch='log',log_a=10000))
    #plt.show()

    ### now make a sample PSF (these are only used for parametric fitting, and we can use WebbPSF to generate one if needed)

    kernel = Gaussian2DKernel(2)
    kernel.normalize()
    psf = kernel.array
    #plt.imshow(psf,origin='lower',cmap='gray')
    #plt.show()


    ### now combine the image with the PSF
    image = convolve(image, psf)
    #plt.imshow(image, origin='lower',norm=simple_norm(image,stretch='log',log_a=10000))
    #plt.show()

    ### applying noise
    np.random.seed(3)
    gain = 1e5
    image = np.random.poisson(image*gain)/gain
    snp = 100.0
    sky_sigma = 1.0/snp
    image+=sky_sigma * np.random.standard_normal(size=(ny,nx))
    #plt.imshow(image,origin='lower',cmap='gray',norm=simple_norm(image,stretch='log',log_a=10000))
    #plt.show()


    #########################################



    ##########################
    """
    Creating the segmentation map from the image
    Use libraries like SExtractor or photutils(use this one)
    """
    ##########################

    # count anything >1.5 sigma as a detection
    threshold = detect_threshold(image, 1.5)
    npixels = 5 # min number of connected pixels
    convolved_image = convolve(image,psf)
    segmap = detect_sources(convolved_image,threshold,npixels)
    #plt.imshow(segmap,origin='lower',cmap='gray')
    #plt.show()
    '''
    ## here value of 0 is the background, and 1 is for the labeled source
    # but statmorph can take in a variety of values for this

    ##########################

    """
    Actually running statmorph

    """
    ##########################

    '''
    start = time.time()
    source_morphs = statmorph.source_morphology(image,segmap, gain=5.8, psf=psf)
    print(f'Time: {time.time()-start}')
    morph = source_morphs[0]


    #### Getting properties from statmorph

    print('Basic measurements (non-parametric)')
    print(f'concentration \t(C): {morph.concentration} ')
    print(f'asymmetry \t(A): {morph.asymmetry} ')
    print(f'smoothness\t(S): {morph.smoothness} ')
    print(f'Gini/M20: {morph.gini}/{morph.m20}')
    print(f'S/N (per pixel): {morph.sn_per_pixel}')


    print('Parametric measurements (Sersic Model)')
    print(f'Sersic amplitude: {morph.sersic_amplitude}')
    print(f'Sersic rhalf: {morph.sersic_rhalf}')
    print(f'Sersic n: {morph.sersic_n}')


    #####################
    """
    Visualizing the fit
    """
    #####################
    #from statmorph.utils.image_diagnostics import make_figure

    #fig = make_figure(morph)

    ### this not working?? why?
    #fig.savefig('tutorial.png')
    #plt.close(fig)

    '''


    #####################################
    """
    Importing a JADES fits file for testing 
    """
    #####################################

    ### It looks like the CAS, Gini/M20 values are the same from statmorph
    ### regardless of gain and psf values
    ### only values to change under this change are SN and parameteric values like sersic


    ### importing fits files for GS ##############

    hdul=None
    gs_whtmap = None
    gs_jades_seg = None
    full_jades_seg = None # for testing with grizli so that we use thE JADES SEGMAP
    gs_jades_wcs = None #   ^^^^^^^^^^^
    ### try to use the grizli data reduction otherwise use jades
    if(use_grizli):
        grizli_hdul = fits.open(f'/home/robbler/research/grizli-reductions/gds-grizli-v7.2-{filter}-clear_drc_sci.fits')
        hdul=grizli_hdul
        gs_whtmap = fits.open(f'/home/robbler/research/grizli-reductions/gds-grizli-v7.2-{filter}-clear_drc_wht.fits')[0].data
        # gs_jades_seg = fits.getdata('/home/robbler/research/grizli-reductions/gds-grizli-v7.2-ir_seg.fits')
        gs_jades_seg = fits.getdata('/home/robbler/research/JADES_Fits_Maps/GS/hlsp_jades_jwst_nircam_goods-s-deep_segmentation_v2.0_drz.fits')
        full_jades_seg = fits.open('/home/robbler/research/JADES_Fits_Maps/GS/hlsp_jades_jwst_nircam_goods-s-deep_segmentation_v2.0_drz.fits')
    else:
        hdul = fits.open(f'/home/robbler/research/JADES_Fits_Maps/GS/hlsp_jades_jwst_nircam_goods-s-deep_{filter}_v2.0_drz.fits')
        # segmentation map
        gs_jades_seg = fits.getdata('/home/robbler/research/JADES_Fits_Maps/GS/hlsp_jades_jwst_nircam_goods-s-deep_segmentation_v2.0_drz.fits')



    # catalog file
    tmp = fits.open('/home/robbler/research/JADES catalog/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits')

    # generate PSF for current filter
    # nc = webbpsf.NIRCam()
    # nc.filter = 'F444W'
    # psf = nc.calc_psf(oversample=4,display=True)
    # print(psf)
    # exit()

    # .fits file (image version)
    if(use_grizli):
        gs_sci = hdul[0].data
        gs_wcs_coords = WCS(hdul[0].header)
        ###### for testing with grizli sources so that we use the JADES SEGMAP instead
        gs_jades_wcs = WCS(full_jades_seg[0].header)
        # gn_wcs_coords = WCS(hdul[0].header)
        full_jades_seg.close()
    else:
        gs_sci = hdul[1].data # for getting the science image
        gs_wcs_coords = WCS(hdul[1].header) # getting wcs from header
    hdul.close()

    # catalog file
    gs_jades_catalog = Table(tmp['FLAG'].data) # extension #2 "FLAG" data
    tmp.close()

    # print(gn_whtmap[0].header['FLAGS_WEIGHT'])
    # exit()

    #################### open goods north files ... ###################################################################################
    # starting with f444 because it should have the clearest source (maybe not clear detail though)
    #hdul = fits.open('JADES_Fits_Maps/hlsp_jades_jwst_nircam_goods-s-deep_f444w_v2.0_drz.fits')
    
    ### try to use the grizli data reduction otherwise use jades
    gn_whtmap = None
    gn_jades_seg = None
    gn_jades_wcs = None
    if(use_grizli):
        grizli_hdul = fits.open(f'/home/robbler/research/grizli-reductions/gdn-grizli-v7.3-{filter}-clear_drc_sci.fits')
        hdul=grizli_hdul
        gn_whtmap = fits.open(f'/home/robbler/research/grizli-reductions/gdn-grizli-v7.3-{filter}-clear_drc_wht.fits')[0].data
        # gn_jades_seg = fits.getdata('/home/robbler/research/grizli-reductions/gdn-grizli-v7.3-ir_seg.fits')
        gn_jades_seg = fits.getdata('research/JADES_Fits_Maps/GN/hlsp_jades_jwst_nircam_goods-n_segmentation_v1.0_drz.fits')
        full_jades_seg = fits.open('research/JADES_Fits_Maps/GN/hlsp_jades_jwst_nircam_goods-n_segmentation_v1.0_drz.fits')
    else:
        hdul = fits.open(f'research/JADES_Fits_Maps/GN/hlsp_jades_jwst_nircam_goods-n_{filter}_v1.0_drz.fits')
        # segmentation map
        gn_jades_seg = fits.getdata('research/JADES_Fits_Maps/GN/hlsp_jades_jwst_nircam_goods-n_segmentation_v1.0_drz.fits')

    # catalog file
    tmp = fits.open('research/JADES catalog/hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits')

    # generate PSF for current filter
    # nc = webbpsf.NIRCam()
    # nc.filter = 'F444W'
    # psf = nc.calc_psf(oversample=4,display=True)
    # print(psf)
    # exit()

    # .fits file (image version)
    if(use_grizli):
        gn_sci = hdul[0].data
        gn_wcs_coords = WCS(hdul[0].header)
        ###### for testing with grizli sources so that we use the JADES SEGMAP instead
        gn_jades_wcs = WCS(full_jades_seg[0].header)
        # gn_wcs_coords = WCS(hdul[0].header)
        full_jades_seg.close()
    else:
        gn_sci = hdul[1].data # for getting the science image
        gn_wcs_coords = WCS(hdul[1].header) # getting wcs from header
    hdul.close()
    # catalog file
    gn_jades_catalog = Table(tmp['FLAG'].data) # extension #2 "FLAG" data
    tmp.close()


    # maybe this is the weight value we're looking for for statmorph (one value for each entry, so how would I convert this to an image to be used with statmorph though?)
    #print(jades_catalog[3]['F444W_WHT'])

    fwhms = {'F070W':0.742,'F090W':0.968,'F115W':1.194,'F150W':1.581,'F200W':2.065,'F277W':1.397,'F356W':1.810,'F444W':2.222}
    # ^ as described in https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-point-spread-functions#gsc.tab=0

    ###################
    """
    Now take the cutout images for all sources
    """
    ###################
    

    sources = get_sources()
    
    # source_names, jades_names, source_coords = get_GS()
    # source_names = np.array(source_names)
        
    # ind = np.where(source_names==kirk_id)[0][0]
    # detection radius ~ 1" for Kirkpatrick data, so we want anything within a ~3" box
    # maybe change this later to reflect the light distribution (20%, 80% radii)

    if(use_grizli):
        size = 8 # cutout size in arcsec, might need to tweak depending on the source size
    # this was a size of 6 for JADES sources (scaling ratio is ~.75)
    else:
        size = 6 # ~0.75 the grizli one (equates to 200x200) sized cutouts

    
    # jades_id = jades_names[ind]
    # get jades id counterpart coordinates
    '''
    #source = jades_catalog[jades_catalog['ID']==jades_id]
    #print(source['ID'])
    position = SkyCoord(source_coords[ind,0],source_coords[ind,1],unit='deg')
    #position = SkyCoord(source['RA'],source['DEC'],unit='deg')
    img = Cutout2D(sci,position,size*units.arcsec,wcs=wcs_coords).data
    source_seg = Cutout2D(jades_seg,position,size*units.arcsec,wcs=wcs_coords).data
    #source_gain = source['F444W_WHT'] # this seems wrong, need to figure out weightmap issue
    source_gain = 100.
    morph = statmorph.source_morphology(img,source_seg,gain=source_gain)[0]

    ################# ^^^^^^^^^^^^^^^^ ###################
    ## the value of gain=source_gain might not be right, need to check with Alex or Sam to see if that's what I should be doing here...


    print(f'Values for JADES:{jades_id}')
    print(f'concentration \t(C): {morph.concentration} ')
    print(f'asymmetry \t(A): {morph.asymmetry} ')
    print(f'smoothness(clumpiness) \t(S): {morph.smoothness} ')
    print(f'Gini: {morph.gini}')
    print(f'M20: {morph.m20}')
    print(f'S/N (per pixel): {morph.sn_per_pixel}')
    print(f'Flag: {morph.flag}')


    # then maybe show the source to see
    fig = make_figure(morph)
    plt.show()
    '''    

    # importing psf from testing fits file first (f444w)
    # psf = fits.getdata(f'/home/robbler/research/JADES catalog/PSF_NIRCam_in_flight_opd_filter_{filter.upper()}.fits')
    ### for testing with a specific ID
    # search_one_id = np.array([kirk_id])
    search_ids = sources['ID']
    full_output = sources
    # search_ids = ['GN_IRS33^e','GS_IRS12','GS_IRS34','GS_IRS37','GS_IRS62','GS_IRS9'] # trouble sources with grizli
    # search_ids = ['GS_IRS2'] # for testing with a high S/N source
    # search_ids = ['GS_IRS34']

    cols = ['Filter','Flag','Concentration (C)','Asymmetry (A)','Outer Asymmetry (Ao)','Smoothness (S)','Gini','M20','Gini-M20-Merger','Multimode (M)','Intensity (I)','Deviation (D)','S/N','Cutout Size (pix)']
    

    ###############################
    #### searching for the ID #####
    ###############################
    for id in search_ids:
        row = []
        if(id=='GS_IRS34'):
            size = 10 # change the size of the cutout for the biggest source
        else:
            size = 8
            # this way it's properly fitting it and not returning flag 2

        data_source = sources[(sources['ID'].str.contains(id,regex=False))].iloc[0]
        search_id = id
        jades_catalog,sci,jades_seg,wcs_coords,jades_wcs_coords,whtmap = (None for i in range(6))

        if(id.startswith('GN')):
            jades_catalog = gn_jades_catalog
            sci = gn_sci
            jades_seg = gn_jades_seg
            wcs_coords = gn_wcs_coords
            jades_wcs_coords = gn_jades_wcs
            whtmap = gn_whtmap
        if(id.startswith('GS')):
            jades_catalog = gs_jades_catalog
            sci = gs_sci
            jades_seg = gs_jades_seg
            wcs_coords = gs_wcs_coords
            jades_wcs_coords = gs_jades_wcs
            whtmap = gs_whtmap

        jades_id_raw = data_source['JADES ID']
        jades_id = 0
        if(' ' in jades_id_raw):
            jades_id = int(jades_id_raw.split(' ')[0]) # splitting because sometimes I have a secondary source on there, but the main is always listed first (this only happens 2 times in our sample)
        else:
            jades_id = int(jades_id_raw)

        source = jades_catalog[jades_catalog['ID']==jades_id] # make sure jades_id is an integer, otherwise it just returns the first item NOT AN ERROR??
        position = SkyCoord(source['RA'],source['DEC'],unit='deg')


        ###################
        # try to get the cutout, print error and continue if it fails #
        # usually fails if the wcs coordinates aren't set right
        ###################
        try:
            img = Cutout2D(sci,position,size*units.arcsec,wcs=wcs_coords,copy=True).data
            wht=None
            if(use_grizli):
                wht = Cutout2D(whtmap,position,size*units.arcsec,wcs=wcs_coords,copy=True).data
                #### testing with different size to use the jades segmap instead so we can get accurate measurement for S/N
                seg = Cutout2D(jades_seg,position,.75*size*units.arcsec,wcs=jades_wcs_coords,copy=True).data
            else:
                # get the segmap for the specific source
                seg = Cutout2D(jades_seg,position,size*units.arcsec,wcs=wcs_coords,copy=True).data
            seg_img = SegmentationImage(seg)
            # plt.imshow(seg_img)
            # plt.show()
            # plt.imshow(img)
            # plt.show()
        except:
            print(f'[!!!] ERROR with IMG for {search_id}')
            continue
############ TESTING SNR #################################
        # sn_testing = img/(1/wht)

        ## testing with something of a similar size to the gini segmap
        # gini_cut = Cutout2D(sci,position,size/7*units.arcsec,wcs=wcs_coords,copy=True).data
        # gini_wht = Cutout2D(whtmap,position,size/7*units.arcsec,wcs=wcs_coords,copy=True).data
        # new_sn = gini_cut/(1/gini_wht)
        # plt.imshow(new_sn,origin='lower')
        # plt.colorbar()
        # plt.show()
        # # print(np.mean(sn_testing))
        # print(np.mean(new_sn))

        # full_sn = img/(1/wht)
        # plt.imshow(full_sn,origin='lower')
        # plt.colorbar()
        # plt.show()
        # print(np.mean(full_sn))
        # exit()
        # # noise = np.std(1/wht)

        # print(img.shape)
        # noise_cut = img[100:150,0:50]
        # plt.imshow(noise_cut,origin='lower')
        # plt.colorbar()
        # plt.show()

        # plt.imshow(img,origin='lower')
        # plt.colorbar()
        # plt.show()
        # noise_new = np.std(noise_cut)

        # print(np.max(img)/noise_new)
        # sn_test = img/noise
        # print(np.max(img)/noise)
        # # print(sn_test)
        # # print(sn_testing)
        # plt.imshow(sn_test,origin='lower')
        # plt.show()
        # exit()
#############################################
        # seg_img_id = np.where(seg_img.labels==jades_id)[0][0] # for running statmorph on only the jades ID we want
        # m = np.mean(img)
        # s = np.std(img)
        # sigma = 2.0 # detection threshold =1.5*sigma
        # fwhm = fwhms[filter.upper()]

        # making a 2d gaussian kernel with fwhm 3 pixels (Ren paper uses 0.2 Petrosian Radius)
        # then getting a sample of the background and using that as the threshold for detection 
        # instead of just sigma ##################
        # npix = 60 # min pixels for connection
        # from astropy.stats import SigmaClip
        # sigma_clip = SigmaClip(sigma=4.0,maxiters=10)
        # sigma_clip = SigmaClip(sigma=sigma+2.0)
        # bkg_estimator = MedianBackground(sigma_clip=sigma_clip)
        # bkg = Background2D(img, (9,9), filter_size=(13,13), bkg_estimator=bkg_estimator)
        # bkg = Background2D(img,(25,25),filter_size=(9,9),sigma_clip=sigma_clip,bkg_estimator=bkg_estimator)
        # img -= bkg.background  # subtract the background
        # kernel = make_2dgaussian_kernel(fwhm, size=5)
        # convolved_data = convolve(img, kernel)
        # convolved_image = convolve(img,psf[0:100][0:100])
        # convolved_image = convolve(img,psf)
        #threshold = detect_threshold(convolved_data,sigma)
        # threshold = detect_threshold(convolved_data,sigma,sigma_clip=sigma_clip)
        #mean,_,std = sigma_clipped_stats(img)
        #threshold = sigma*std
        # make the segmentation map using the bg subtracted image, threshold, & number of required connected pixels
        # segmap = detect_sources(convolved_data,threshold,npix)
        # then separate connected segmentation sources by deblending them
        # deblended_segmap = deblend_sources(convolved_data,segmap,npixels=npix,nlevels=32,contrast=0.5) # nlevel=32 & contrast=0.01
        
        # segmap_final = convolve(deblended_segmap,seg)
        # deblended_segmap = seg # testing with no deblending
        # plt.imshow(segmap,origin='lower')
        # plt.show()
        # from matplotlib.gridspec import GridSpec
        # fig = plt.figure(figsize=(10,6))
        # gs = GridSpec(1, 2, figure=fig)
        # fig.add_subplot(gs[0, 0]).imshow(deblended_segmap,origin='lower')
        # fig.add_subplot(gs[0, 1]).imshow(img,origin='lower',cmap='gray_r')
        # plt.show()
        # exit()
        # plt.imshow(img, origin='lower',norm=colors.PowerNorm(gamma=0.5,vmin=(m-s),vmax=(m+8*s)),cmap='gray_r')
        # plt.show()

        ## attempt to find the biggest segmap part, and fit that within the cutout
        # areas = np.zeros([len(seg_img.labels)])
        # v=0
        # areas = seg_img.get_areas(seg_img.labels)
        # # for i in seg_img.labels:
        #     # v+=1
        # seg_id = seg_img.labels[np.where(areas==np.max(areas))[0][0]]
        #################################################

        source_gain = source[f'{filter.upper()}_WHT'] # this seems wrong, need to figure out weightmap issue
        # source_gain = 100
        # source_gain = 1.82*60000 #(source: https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-detector-overview/nircam-detector-performance#gsc.tab=0)
        # source_gain = 1e9
        ### ^^^^ maybe useful later for all nircam data gain map
        ### https://jwst-pipeline.readthedocs.io/en/latest/jwst/gain_scale/description.html
        ## the value of gain=source_gain might not be right, need to check with Alex or Sam to see if that's what I should be doing here...
        print(f'----------------')
        morph = None
        if(use_grizli):
            morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=(1/wht),skybox_size=32,verbose=True) # weight map some auto generated by SExtractor will give it as 1/RMS
            # maybe label id should be seg_id (but it's only getting the max value here)
        else:
            morph = statmorph.SourceMorphology(img,seg_img,jades_id,gain=source_gain,skybox_size=32,verbose=True) # maybe cutout_extent=1?
        
        print(f'[---] Searching for Kirk ID: {search_id}')
        print(f'[---] Searching for JADES ID: {jades_id}')
        print(f'[---] Found JADES ID: {source["ID"][0]}')
        print(f'[>>>] Cutout size being inputted: {img.shape}')
        print(f'Source location: {position}')
        
        row.append(filter)
        row.append(morph.flag)
        row.append(morph.concentration)
        row.append(morph.asymmetry)
        row.append(morph.outer_asymmetry)
        row.append(morph.smoothness)
        row.append(morph.gini)
        row.append(morph.m20)
        row.append(morph.gini_m20_merger)
        row.append(morph.multimode)
        row.append(morph.intensity)
        row.append(morph.deviation)
        row.append(morph.sn_per_pixel)
        row.append(str(img.shape))

        full_output.loc[full_output['ID']== search_id, cols] = row

        print(f'Flag: {morph.flag}')
        print(f'concentration \t(C): {morph.concentration} ')
        print(f'asymmetry (under-estimate) (A): {morph.asymmetry} ')
        print(f'outer asymmetry\t (A0): {morph.outer_asymmetry}')
        print(f'smoothness(depends on psf heavily) (S): {morph.smoothness} ')
        print(f'Gini: {morph.gini}')
        print(f'M20: {morph.m20}')
        print(f'Gini-M20 merger stat: {morph.gini_m20_merger}')
        print(f'Multimode (M): {morph.multimode}')
        print(f'Intensity: {morph.intensity}')
        print(f'Deviation: {morph.deviation}')
        print(f'S/N (per pixel): {morph.sn_per_pixel}')
        
        fig = make_figure(morph)
        savepath = (f'research/statmorph_output/{filter}/{size}_{search_id}')
        if(use_grizli):
            savepath = (f'research/statmorph_output/grizli/{filter}/{size}_{search_id}')
        plt.savefig(savepath,dpi=100)
        # plt.show()
        plt.close()

    full_output.to_csv(f'research/statmorph_output/grizli/{filter}_{size}_{output_name}.tsv','\t')
    print('[!!!] Finished statmorph measurements -- outputing table...')
    print(full_output)
    disclaimer()
    return

##################
"""
::: NEXT STEPS ::::

Confirm this is working properly by looking at the values from other data(?)
Confirm we're getting all the values we want from this (ensure that we don't want reff or n or other parametric values)
Then run this for all sources in the sample
and find a way to save all this to a pandas database or maybe a astropy table
"""
##################

def test_in_ceers():
    '''test it in the ceers survey to reference results with this paper: https://arxiv.org/pdf/2404.16686'''
    #catalog_filename = '/home/robbler/research/CEERS_statmorph_testing/ceers5_f200w_cat.ecsv'
    catalog_filename = '/home/robbler/research/CEERS_statmorph_testing/hlsp_candels_hst_wfc3_egs_multi_v1_mass-cat.fits'
    table = Table()
    with fits.open(catalog_filename) as hdul:
        # this way it closes automatically when with: ends
        hdr = hdul[0].header
        flag_data = Table(hdul[1].data)
        table = flag_data['ID','RAdeg','DECdeg']
        search_ids = np.array([13447,21320,19427,6944,6611,31265,19511,13206,2705,14053,20998,918,31485,3747,11689,20248,11790,23584,23419,13219,16315,12403,2985])
        #search_id = 21320
        c=0

        ###############################################
        #### creating a table just to get practice ####
        ###############################################
        rows = np.zeros((len(search_ids),3))


        for i in search_ids:
            mask = table['ID']==i
            rows[c,1] = table[mask]['RAdeg']
            rows[c,0] = i
            rows[c,2] = table[mask]['DECdeg']
            c+=1


        table = Table(rows=rows, names=['id', 'ra','dec'])
        #make_circles(table,1.0)
        #exit()
        # data['dec'] = 
        
    

    #### now import the fits file and seg map for statmorph

    #hdul = fits.open('JADES_Fits_Maps/hlsp_jades_jwst_nircam_goods-s-deep_f444w_v2.0_drz.fits')
    #search_id = 6611

    search_ids = np.array([23419,21320,31265,31485,23584,13447,19511,13219,16315,11689,12403,11790,13206,19427,20248,20998,14053,6944,3747,2705,2985,6611,918])
    fits_loc = np.array([1,2,2,3,3,4,4,4,4,4,5,6,6,6,6,6,6,8,8,8,10,10,7])
    fits_dict = {
        1: "/home/robbler/research/CEERS_statmorph_testing/hlsp_ceers_jwst_nircam_nircam1_f200w_v0.5_i2d.fits.gz",
        2: "/home/robbler/research/CEERS_statmorph_testing/hlsp_ceers_jwst_nircam_nircam2_f200w_v0.5_i2d.fits.gz",
        3: "/home/robbler/research/CEERS_statmorph_testing/hlsp_ceers_jwst_nircam_nircam3_f200w_v0.5_i2d.fits.gz",
        4: "/home/robbler/research/CEERS_statmorph_testing/hlsp_ceers_jwst_nircam_nircam4_f200w_dr0.6_i2d.fits.gz",
        5: "/home/robbler/research/CEERS_statmorph_testing/hlsp_ceers_jwst_nircam_nircam5_f200w_dr0.6_i2d.fits.gz",
        6: "/home/robbler/research/CEERS_statmorph_testing/hlsp_ceers_jwst_nircam_nircam6_f200w_v0.5_i2d.fits.gz",
        7: "/home/robbler/research/CEERS_statmorph_testing/hlsp_ceers_jwst_nircam_nircam7_f200w_dr0.6_i2d.fits.gz",
        8: "/home/robbler/research/CEERS_statmorph_testing/hlsp_ceers_jwst_nircam_nircam8_f200w_dr0.6_i2d.fits.gz",
        9: "/home/robbler/research/CEERS_statmorph_testing/hlsp_ceers_jwst_nircam_nircam9_f200w_dr0.6_i2d.fits.gz",
        10: "/home/robbler/research/CEERS_statmorph_testing/hlsp_ceers_jwst_nircam_nircam10_f200w_dr0.6_i2d.fits.gz"
    }

    ##### run loop on all them ######
    statmorph_table = Table()
    stat_rows = np.zeros((len(search_ids),13),dtype=object)
    c=0

    ### for testing with a specific ID
    search_one_id = np.array([11689])
    for id in search_one_id:
    # for id in search_ids:
        search_id = id
        fits_path = fits_dict[fits_loc[np.where(search_ids == id)[0][0]]]
        hdul = fits.open(fits_path)
        #hdul = fits.open('/home/robbler/research/CEERS_statmorph_testing/hlsp_ceers_jwst_nircam_nircam10_f200w_dr0.6_i2d.fits.gz')
        # segmentation map
        # don't have the segmentation map available so need to use something like WebbPSF to find psf then create a segmentation map 
        # from this with anything that has a SNR of >2.5 (standard) or sigma of ??? 
        #seg = fits.getdata('/home/robbler/research/CEERS_statmorph_testing/ceers5_f200w_segm.fits.gz')
        # importing psf
        psf = fits.getdata('/home/robbler/research/CEERS_statmorph_testing/ceers-full-grizli-v6.0-f200w-clear_drc_sci.gz_psf.fits')
        sci = hdul[1].data # for getting the science image
        wcs_coords = WCS(hdul[1].header) # getting wcs from header
        hdul.close()
        size = 3.0 # cutout size in arcsec, might need to tweak depending on the source size
        mask = table['id']==search_id
        print(f'searching for id: {search_id}')
        position = SkyCoord(table[mask]['ra'],table[mask]['dec'],unit='deg')
        img = Cutout2D(sci,position,size*units.arcsec,wcs=wcs_coords,copy=True).data
        #source_seg = Cutout2D(seg,position,size*units.arcsec,wcs=wcs_coords,copy=True).data
        m = np.mean(img)
        s = np.std(img)
        sigma = 1.5
        fwhm = 10.0
        # making a 2d gaussian kernel with fwhm 3 pixels (Ren paper uses 0.2 Petrosian Radius)
        # then getting a sample of the background and using that as the threshold for detection 
        # instead of just sigma ##################
        npix = 5 # min pixels for connection
        #bkg_estimator = MedianBackground()
        # from astropy.stats import SigmaClip
        # sigma_clip = SigmaClip(sigma=3.5)
        #bkg = Background2D(img, (9,9), filter_size=(21,21), bkg_estimator=bkg_estimator)
        #img -= bkg.background  # subtract the background
        #kernel = make_2dgaussian_kernel(fwhm, size=5)
        #convolved_data = convolve(img, kernel)
        # convolved_image = convolve(convolved_data,psf)
        convolved_image = convolve(img,psf)
        #threshold = detect_threshold(convolved_image,sigma)
        from astropy.stats import sigma_clipped_stats
        mean,_,std = sigma_clipped_stats(img)
        threshold = sigma*std
        # make the segmentation map using the bg subtracted image, threshold, & number of required connected pixels
        segmap = detect_sources(convolved_image,threshold,npix)
        
        # then separate connected segmentation sources by deblending them
        deblended_segmap = deblend_sources(convolved_image, segmap,npixels=npix,nlevels=32,contrast=0.01)
        # deblended_segmap = segmap # testing with no deblending
        plt.imshow(deblended_segmap,origin='lower')
        plt.show()
        # plt.imshow(img, origin='lower',norm=colors.PowerNorm(gamma=0.5,vmin=(m-s),vmax=(m+8*s)),cmap='gray_r')
        # plt.show()

        ## attempt to find the biggest segmap part, and fit that within the cutout
        areas = np.zeros([len(deblended_segmap.labels)])
        for i in deblended_segmap.labels:
            areas[i-1] = deblended_segmap.get_area(i)
        seg_id = np.where(areas==np.max(areas))[0][0]
        #################################################

        #source_gain = source['F444W_WHT'] # this seems wrong, need to figure out weightmap issue
        source_gain = 1e5 # used in example (apparently fairly high though)
        ### ^^^^ maybe useful later for all nircam data gain map
        ### https://jwst-pipeline.readthedocs.io/en/latest/jwst/gain_scale/description.html
        print(f'Cutout size: {img.shape}')
        morph_full = statmorph.source_morphology(img,deblended_segmap,gain=source_gain,psf=psf,cutout_extent=1.5,verbose=True) # maybe cutout_extent=1?
        morph = morph_full[seg_id]
        
        ###### ^ so it's fitting multiple of them depending on which label in the segmentation map it is ^^^^^^
        print(f'morph label (from segmap): {morph.label}, total # of morphs:{len(morph_full)}')
        ## https://github.com/vrodgom/statmorph/blob/master/docs/description.rst
        ################# ^^^^^^^^^^^^^^^^ ###################
        ## the value of gain=source_gain might not be right, need to check with Alex or Sam to see if that's what I should be doing here...
        stat_rows[c,0] = id #id
        stat_rows[c,1] = position # position (skycoord)
        stat_rows[c,2] = morph.flag # flag (error, 0=great, 1=ok, 2=bad)
        stat_rows[c,3] = morph.concentration # concentration
        stat_rows[c,4] = morph.asymmetry # asymmetry
        stat_rows[c,5] = morph.outer_asymmetry # outer asymmetry
        stat_rows[c,6] = morph.deviation # deviation (not sure if this is total or just outer?)
        stat_rows[c,7] = morph.smoothness # smoothness (clumpiness)
        stat_rows[c,8] = morph.gini # gini
        stat_rows[c,9] = morph.m20 # m20
        stat_rows[c,10] = morph.gini_m20_merger # gini/m20 for mergers
        stat_rows[c,11] = morph.sn_per_pixel # s/n ratio (per pixel?)
        stat_rows[c,12] = (f'({morph.xmax_stamp-morph.xmin_stamp},{morph.ymax_stamp-morph.ymin_stamp})') # size of cutout
        c+=1

        '''
        print(f'Values for CANDELS(EGS):{search_id}')
        print(f'Source location: {position}')
        print(f'Flag (0=great,1=ok,2=bad): {morph.flag}')
        print(f'concentration \t(C): {morph.concentration} ')
        print(f'asymmetry (under-estimate) (A): {morph.asymmetry} ')
        print(f'outer asymmetry\t (A0): {morph.outer_asymmetry}')
        print(f'Deviation (not outer?) (D): {morph.deviation}')
        print(f'smoothness(depends on psf heavily) (S): {morph.smoothness} ')
        print(f'Gini: {morph.gini}')
        print(f'M20: {morph.m20}')
        print(f'Gini-M20 merger stat: {morph.gini_m20_merger}')
        print(f'S/N (per pixel): {morph.sn_per_pixel}')
        print(f'Size of sample cutout (WCS): ({morph.xmax_stamp-morph.xmin_stamp},{morph.ymax_stamp-morph.ymin_stamp})')
        '''
        
        fig = make_figure(morph)
        plt.show()
        # plt.savefig(f'/home/robbler/research/CEERS_statmorph_testing/statmorph_fits/{id}.png')
        # plt.close()
    statmorph_table = Table(rows=stat_rows, names=['id', 'position','flag','concentration (C)','asymmetry (A)','outer asymmetry (Ao)','deviation','smoothness','gini','m20','gini/m20 merger','s/n ratio','cutout size'])
    table_output = statmorph_table.to_pandas()
    table_output.to_csv("/home/robbler/research/CEERS_statmorph_testing/statmorph_measurements.csv", encoding='utf-8', index=False)



# test_in_ceers()

grizli_filters = ['f115w','f150w','f200w','f277w','f356w','f444w']
grizli_filters = ['f277w','f356w','f444w'] # testing with these bc others are giving errors?
jades_filters = ['f115w','f150w','f200w','f277w','f356w','f444w']
for i in grizli_filters:
    run_in_jades(filter=i,use_grizli=True,output_name='grizli_1_over_wht')





def make_circles(table,radius):
    '''function to make regions file given equal size ra,dec,label arrays with size radius.
        ra & dec can be in (") or (deg)
        radius in arcseconds(")'''
    c = 0
    print("""# Region file format: DS9 version 4.1
global color=cyan dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5""")
    for i in table:
        label = i['id']
        ra = i['ra']
        dec = i['dec']
        print(f'circle({ra}, {dec}, {radius}")  # text={{{label}}}')
        c+=1


# save all output to "all, err, out" files
# command 2> >(tee /research/statmorph_output/err) 1> >(tee /research/statmorph_output/out) | tee >/research/statmorph_output/all.txt