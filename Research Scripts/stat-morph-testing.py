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
from scipy.ndimage import zoom


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


def run_in_jades(output_name='Grizli-statmorph-measurements',use_grizli=False,error_bars=False,source_search_id="",source_in_filter=""):
    ###############################
    """
    Get statmorph fit for a jades source
    Input:  output_name: Filename to output the table of measurements to (saved filename is inclusive of filter name)
            use_grizli: Use Grizli Reductions or not, if false, defaults to using JADES reductions (less data available)
            error_bars: Apply 100 Gaussians to the input image to statmorph to get 16th & 84th percentile from that for error bars
                        on the following measurements: Asymmetry, Concentration, Gini, M20 -- others must be added manually
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
    #####################################
    """
    Importing a JADES fits file for testing 
    """
    #####################################


    ### importing fits files for GS ##############

    ######## store everything in respective array location for the following Nircam filters
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
        if(use_grizli):
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
        else:
            hdul = fits.open(f'/home/robbler/research/JADES_Fits_Maps/GS/hlsp_jades_jwst_nircam_goods-s-deep_{f}_v2.0_drz.fits')
            # segmentation map
            gs_jades_seg = fits.getdata('/home/robbler/research/JADES_Fits_Maps/GS/hlsp_jades_jwst_nircam_goods-s-deep_segmentation_v2.0_drz.fits')
            gs_scis[ind] = hdul[1].data # for getting the science image
            gs_wcs_coords[ind] = WCS(hdul[1].header) # getting wcs from header


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

        if(use_grizli):
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
        else:
            hdul = fits.open(f'research/JADES_Fits_Maps/GN/hlsp_jades_jwst_nircam_goods-n_{f}_v1.0_drz.fits')
            # segmentation map
            gn_jades_seg = fits.getdata('research/JADES_Fits_Maps/GN/hlsp_jades_jwst_nircam_goods-n_segmentation_v1.0_drz.fits')
            gn_scis[ind] = hdul[1].data # for getting the science image
            gn_wcs_coords[ind] = WCS(hdul[1].header) # getting wcs from header

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

    # maybe this is the weight value we're looking for for statmorph (one value for each entry, so how would I convert this to an image to be used with statmorph though?)
    #print(jades_catalog[3]['F444W_WHT'])

    fwhms = {'F070W':0.742,'F090W':0.968,'F115W':1.194,'F150W':1.581,'F200W':2.065,'F277W':1.397,'F356W':1.810,'F444W':2.222}
    # ^ as described in https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-point-spread-functions#gsc.tab=0

    ###################
    """
    Now take the cutout images for all sources
    """
    ###################
    

    ### testing a new whtmap value from https://dawn-cph.github.io/dja/blog/2023/07/18/image-data-products/
    # import astropy.io.fits as pyfits
    # import scipy.ndimage as nd
    # sci_img = pyfits.open('/home/robbler/research/grizli-reductions/gdn-grizli-v7.3-f356w-clear_drc_sci.fits')
    # # wht_img = pyfits.open('/home/robbler/research/grizli-reductions/gdn-grizli-v7.3-f356w-clear_drc_wht.fits')
    # exp_img = pyfits.open('/home/robbler/research/grizli-reductions/gdn-grizli-v7.3-f356w-clear_drc_exp.fits')
    # full_exp = np.zeros(sci_img[0].data.shape,dtype=int)
    # full_exp[2::4,2::4] += exp_img[0].data*1
    # full_exp = nd.maximum_filter(full_exp,4)

    # # full_exp_img = pyfits.HDUList([pyfits.PrimaryHDU(data=full_exp)])[0].data
    # full_exp_img = pyfits.HDUList([pyfits.PrimaryHDU(data=full_exp)])

    # full_exp_img.writeto('/home/robbler/research/grizli-reductions/full_img_exposure_f356w.fits')

    # exit()

    # full_exp_img = fits.open('/home/robbler/research/grizli-reductions/full_img_exposure_f356w.fits')[0].data
    
    sources = get_sources(full=True)
    
    # source_names, jades_names, source_coords = get_GS()
    # source_names = np.array(source_names)
        
    # ind = np.where(source_names==kirk_id)[0][0]
    # detection radius ~ 1" for Kirkpatrick data, so we want anything within a ~3" box
    # maybe change this later to reflect the light distribution (20%, 80% radii)

    # importing psf from testing fits file first (f444w)
    # psf = fits.getdata(f'/home/robbler/research/JADES catalog/PSF_NIRCam_in_flight_opd_filter_{filter.upper()}.fits')
    ### for testing with a specific ID
    # search_one_id = np.array([kirk_id])
    search_ids = sources['ID']
    full_output = sources
    # search_ids = ['GN_IRS33^e','GS_IRS12','GS_IRS34','GS_IRS37','GS_IRS62','GS_IRS9'] # trouble sources with grizli
    # search_ids = ['GS_IRS2'] # for testing with a high S/N source
    # search_ids = ['GN_IRS2','GN_IRS21','GN_IRS22^d','GN_IRS25','GN_IRS33^e','GN_IRS4','GN_IRS55','GN_IRS56','GN_IRS60','GN_IRS61','GS_IRS20','GS_IRS37','GS_IRS50','GS_IRS60']
    # search_ids = ['GS_IRS12','GN_IRS55','GS_IRS73']


    # if(source_search_id!=""):
    #     search_ids = [i for i in source_search_id.split(',')]
    ### search ids for testing error bar distribution
    # non f277w rest-frame
    search_ids = ['GN_IRS14','GN_IRS17','GN_IRS2','GN_IRS24^e','GN_IRS27','GN_IRS33^e','GN_IRS36','GN_IRS4','GN_IRS43','GN_IRS60','GN_IRS61','GN_IRS8^d','GS_IRS2','GS_IRS21','GS_IRS25','GS_IRS37','GS_IRS62','GS_IRS70^e','GS_IRS71^e','GS_IRS72^e','GS_IRS73','GS_IRS9']
    # f277w rest-frame
    # search_ids = ['GN_IRS10','GN_IRS11','GN_IRS15','GN_IRS21','GN_IRS22^d','GN_IRS25','GN_IRS55','GN_IRS56','GN_IRS58^d','GN_IRS59^d','GN_IRS7^d','GS_IRS12','GS_IRS14','GS_IRS15','GS_IRS20','GS_IRS23','GS_IRS34','GS_IRS45','GS_IRS50','GS_IRS58^f','GS_IRS60','GS_IRS61','GS_IRS64^e']
    # search_ids = ['GS_IRS12',
    # 'GS_IRS14',
    # 'GS_IRS15',
    # 'GS_IRS2',
    # 'GS_IRS20',
    # 'GS_IRS21',
    # 'GS_IRS23',
    # 'GS_IRS25',
    # 'GS_IRS34',
    # 'GS_IRS37',
    # 'GS_IRS45',
    # 'GS_IRS50',
    # 'GS_IRS58^f',
    # 'GS_IRS60',
    # 'GS_IRS61',
    # 'GS_IRS62',
    # 'GS_IRS64^e',
    # 'GS_IRS70^e',
    # 'GS_IRS71^e',
    # 'GS_IRS72^e',
    # 'GS_IRS73',
    # 'GS_IRS9']
    

    ### GN sources:
    '''
    'GN_IRS10',
    'GN_IRS11',
    'GN_IRS14',
    'GN_IRS15',
    'GN_IRS17',
    'GN_IRS2',
    'GN_IRS21',
    'GN_IRS22^d',
    'GN_IRS24^e',
    'GN_IRS25',
    'GN_IRS27',
    'GN_IRS33^e',
    'GN_IRS36',
    'GN_IRS4',
    'GN_IRS43',
    'GN_IRS55',
    'GN_IRS56',
    'GN_IRS58^d',
    'GN_IRS59^d',
    'GN_IRS60',
    'GN_IRS61',
    'GN_IRS7^d',
    'GN_IRS8^d'
    '''
    
    ### GS SOURCES: 
    '''
    'GS_IRS12',
    'GS_IRS14',
    'GS_IRS15',
    'GS_IRS2',
    'GS_IRS20',
    'GS_IRS21',
    'GS_IRS23',
    'GS_IRS25',
    'GS_IRS34',
    'GS_IRS37',
    'GS_IRS45',
    'GS_IRS50',
    'GS_IRS58^f',
    'GS_IRS60',
    'GS_IRS61',
    'GS_IRS62',
    'GS_IRS64^e',
    'GS_IRS70^e',
    'GS_IRS71^e',
    'GS_IRS72^e',
    'GS_IRS73',
    'GS_IRS9'
    '''



    
    cols = ['Filter','Flag','Concentration (C)','Concentration (C) Error (16%)','Concentration (C) Error (84%)','Concentration (C) Error Full','Asymmetry (A)','Asymmetry (A) Error (16%)','Asymmetry (A) Error (84%)','Asymmetry (A) Error Full','Outer Asymmetry (Ao)','Smoothness (S)','Gini','Gini Error (16%)','Gini Error (84%)','Gini Error Full','M20','M20 Error (16%)','M20 Error (84%)','M20 Error Full','Gini-M20-Merger','Multimode (M)','Intensity (I)','Deviation (D)','S/N','Cutout Size (pix)']
    ## size dictionary for the arcsec cutout for each
    ## some sources are smaller some are larger etc
    size_dict = {
'GN_IRS10':8,
'GN_IRS11':8,
'GN_IRS14':6,
'GN_IRS15':8,
'GN_IRS17':6,
'GN_IRS2':5,
'GN_IRS21':6,
'GN_IRS22^d':5,
'GN_IRS24^e':5,
'GN_IRS25':8,
'GN_IRS27':8,
'GN_IRS33^e':5,
'GN_IRS36':8, # maybe 7 instead?
'GN_IRS4':5,
'GN_IRS43':5,
'GN_IRS55':6,
'GN_IRS56':7,
'GN_IRS58^d':7,
'GN_IRS59^d':7,
'GN_IRS60':5,
'GN_IRS61':5,
'GN_IRS7^d':7,
'GN_IRS8^d':8,
'GS_IRS12':8,
'GS_IRS14':8,
'GS_IRS15':8,
'GS_IRS2':5,
'GS_IRS20':7,
'GS_IRS21':5,
'GS_IRS23':8,
'GS_IRS25':5,
'GS_IRS34':10,
'GS_IRS37':8, #source not in grizli reductions
'GS_IRS45':8,
'GS_IRS50':7,
'GS_IRS58^f':7, # trying 7, 8 works fine
'GS_IRS60':7,
'GS_IRS61':7,
'GS_IRS62':6,
'GS_IRS64^e':6,
'GS_IRS70^e':5, #trying 5, 6 works fine
'GS_IRS71^e':6,
'GS_IRS72^e':5,
'GS_IRS73':6,
'GS_IRS9':7,
}
    #### search for a specific filter to make sure it runs that filter right

    ###############################
    #### searching for the ID #####
    ###############################
    id_num = 0
    for id in search_ids:
        row = []
        try:
            data_source = sources[(sources['ID'].str.contains(id,regex=False))].iloc[0]
        except:
            print(f'[!!!] {id} not found in data_import sheet, check import_data.py...')
            continue
        source_row = sources.loc[sources['ID'] == id]
        source_filter = source_row['Obs Filter'].iloc[0]
        if(source_in_filter!=""):
            source_filter = source_in_filter
        index = nir_filters.index(source_filter)

        # if(id=='GS_IRS34' and use_grizli): 
        #     size = 10 # change the size of the cutout for the biggest source
        # elif(use_grizli and (source_filter=='f150w' or source_filter=='f200w')):
        #     ## set the size to 1/2 the usual because scaling is ~2x bigger for these filters
        #     ## for GS_irs34 above, it doesn't matter bc its f277w
        #     size = 6
        # elif(use_grizli):
        #     size = 8
        #     # cutout size in arcsec, might need to tweak depending on the source size
        #     # this was a size of 6 for JADES sources (scaling ratio is ~.75)
        # else: # if we're using JADES data only
        #     size = 6 # ~0.75 the grizli one (equates to 200x200) sized cutouts 
        

        ### use the size testing values from size_dict, individually selected sizes for different sized galaxies
        size = size_dict[id]


        search_id = id
        jades_catalog,sci,jades_seg,wcs_coords,jades_wcs_coords,whtmap = (None for i in range(6))
        if(id.startswith('GN')):
            jades_catalog = gn_jades_catalog
            jades_sizes = jades_gn_sizes
            sci = gn_scis[index]
            jades_seg = gn_jades_seg
            wcs_coords = gn_wcs_coords[index]
            exp_wcs_coord = gn_exp_wcs_coords[index]
            wht_wcs = gn_wht_wcs[index]
            jades_wcs_coords = gn_jades_wcs
            whtmap = gn_whtmaps[index]
            expmap = gn_expmaps[index]

            #for poisson noise
            PHOTMJSR = gn_noise_vars[index][0]
            PHOTSCAL = gn_noise_vars[index][1]
            PHOTFNU = gn_noise_vars[index][2]
            OPHOTFNU = gn_noise_vars[index][3]
        if(id.startswith('GS')):
            jades_catalog = gs_jades_catalog
            jades_sizes = jades_gs_sizes
            sci = gs_scis[index]
            jades_seg = gs_jades_seg
            wcs_coords = gs_wcs_coords[index]
            exp_wcs_coord = gs_exp_wcs_coords[index]
            wht_wcs = gs_wht_wcs[index]
            jades_wcs_coords = gs_jades_wcs
            whtmap = gs_whtmaps[index]
            expmap = gs_expmaps[index]
            PHOTMJSR = gs_noise_vars[index][0]
            PHOTSCAL = gs_noise_vars[index][1]
            PHOTFNU = gs_noise_vars[index][2]
            OPHOTFNU = gs_noise_vars[index][3]

        jades_id_raw = data_source['JADES ID']
        jades_id = 0
        if(' ' in jades_id_raw):
            jades_id = int(jades_id_raw.split(' ')[0]) # splitting because sometimes I have a secondary source on there, but the main is always listed first (this only happens 2 times in our sample)
        else:
            jades_id = int(jades_id_raw)

        source = jades_catalog[jades_catalog['ID']==jades_id]
        # source_size = jades_sizes[jades_sizes['ID']==jades_id]
        position = SkyCoord(source['RA'],source['DEC'],unit='deg')
        # ### set the size to something comparable to the segmap size?
        
        # seg_img = SegmentationImage(Cutout2D(jades_seg,position,1.5*10*units.arcsec,wcs=jades_wcs_coords,copy=True).data)
        # print(seg_img.get_area(jades_id)/5000)
        # size = int(seg_img.get_area(jades_id)/5000)
        # # exit()
        ###################
        # try to get the cutout, print error and continue if it fails #
        # usually fails if the wcs coordinates aren't set right
        ###################
        try:
            img = Cutout2D(sci,position,size*units.arcsec,wcs=wcs_coords,copy=True).data
            wht=None
                # full_exp_cutout = None
            if(use_grizli):
                wht = Cutout2D(whtmap,position,size*units.arcsec,wcs=wht_wcs,copy=True).data
                ### different scaling for GS exp map
                ### but this visually looks off when I plot it, like it looks too big of a region
                ### the only scale in the header is the 'sample rate'=4 so maybe need to just resize the image
                ### maybe use opencv and/or skikit image to resize it with interpolation?
                ### https://www.quora.com/How-do-I-resize-interpolate-a-numpy-array-with-floating-numbers
                # if(id.startswith('GS')):
                    # exp = Cutout2D(expmap,position,size*4*units.arcsec,wcs=exp_wcs_coord,copy=True).data
                # else:
                exp = Cutout2D(expmap,position,size*units.arcsec,wcs=exp_wcs_coord,copy=True).data
                seg = Cutout2D(jades_seg,position,size*units.arcsec,wcs=jades_wcs_coords,copy=True).data    
                # resize the segmap to the image size to get proper pixel scaling
                new_seg_size = img.shape
                new_seg = zoom(seg, (new_seg_size[0] / seg.shape[0], new_seg_size[1] / seg.shape[1]), order=0) 
                seg=new_seg
                # if(source_filter=='f150w' or source_filter=='f200w'):
                #     # different scaling for segmentation map on 150 & 200
                #     seg = Cutout2D(jades_seg,position,1.5*size*units.arcsec,wcs=jades_wcs_coords,copy=True).data
                # else:
                #     seg = Cutout2D(jades_seg,position,0.75*size*units.arcsec,wcs=jades_wcs_coords,copy=True).data
                # # seg = Cutout2D(jades_seg,position,0.75*size*units.arcsec,wcs=jades_wcs_coords,copy=True).data

                # #### testing full_exp for weight map instead because it includes the source poisson noise?
                # full_exp_cutout = Cutout2D(full_exp_img,position,size*units.arcsec,wcs=wcs_coords,copy=True).data
            else:
                # get the segmap for the specific source
                seg = Cutout2D(jades_seg,position,size*units.arcsec,wcs=wcs_coords,copy=True).data
            seg_img = SegmentationImage(seg)
                    # plt.imshow(seg_img)
                    # plt.show()
                    # exit()
                    # plt.imshow(img)
                    # plt.show()
                    # plt.imshow(wht)
                    # plt.show()
                    # exit()
        except:
            print(f'[!!!] ERROR with IMG for {search_id}')
            continue


############ TESTING SNR #################################




        #### testing with the masking feature...

        # sn_testing = img/(1/full_exp_cutout)

        ## testing with something of a similar size to the gini segmap

        # sn_size = size/9.33*units.arcsec
        # noise_pos = SkyCoord(189.2365671,62.1340002,unit='deg')
        # pure_noise = Cutout2D(sci/whtmap,noise_pos,sn_size,wcs=wcs_coords,copy=True).data
        # noise_sum = np.sum(pure_noise)
        
        
        # gini_cut = Cutout2D(sci,position,sn_size,wcs=wcs_coords,copy=True).data

        # # print(f'S/N: {np.sum(gini_cut)/noise_sum}')

        # gini_wht = Cutout2D(whtmap,position,sn_size,wcs=wcs_coords,copy=True).data

        # gini_exp = Cutout2D(full_exp_img,position,sn_size,wcs=wcs_coords,copy=True).data
        # new_sn = gini_cut/(1/gini_exp)
        # plt.imshow(gini_cut,origin='lower')
        # plt.colorbar()
        # plt.show()
        # plt.imshow(gini_wht,origin='lower')
        # plt.colorbar()
        # plt.show()
        # exit()


        if(0. in img or 0. in wht): ### checking for 0 values in img, if any change them to np.nan (seems to be issue only with grizli reductions)
            ### something about masking the effects of stars, but they're not stars so its about it being automated??
            print('---- errrrrrrr 0 in img or wht for source -----')
            info = source_row[['ID','Obs Filter']]
            print(f'{info}')

            img[np.where(img==0)] = np.nan
            wht[np.where(wht==0)] = np.nan


        # continue
        # wht_sn = gini_cut/np.sqrt(1/gini_wht)
        # plt.imshow(wht_sn,origin='lower')
        # plt.colorbar()
        # plt.show()
        # exit()

        # plt.imshow(img/(1/wht),origin='lower')
        # plt.colorbar()
        # plt.show()
        # exit()
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

        source_gain = source[f'{source_filter.upper()}_WHT'] # this seems wrong, need to figure out weightmap issue
        # source_gain = 100
        # source_gain = 1.82*60000 #(source: https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-detector-overview/nircam-detector-performance#gsc.tab=0)
        # source_gain = 1e9
        ### ^^^^ maybe useful later for all nircam data gain map
        ### https://jwst-pipeline.readthedocs.io/en/latest/jwst/gain_scale/description.html
        ## the value of gain=source_gain might not be right, need to check with Alex or Sam to see if that's what I should be doing here...
        print(f'----------------')
        morph = None

        
        ###############################################################################
        ################# testing with Gaussian Noise for error bars ##################
        # root_morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=np.sqrt(1/wht),skybox_size=32,verbose=True)
        # fig = make_figure(root_morph)
        # plt.show()
        # plt.close()
        # # snp = 100.0 # seems to work the best, but not sure what it does necessarily
        # snp = root_morph.sn_per_pixel
        # sky_sigma = 1.0 / snp
        # morphs = np.zeros(100,dtype=object)
        # for i in range(100):
        #     img_noise = img+sky_sigma * np.random.standard_normal(size=(img.shape[0],img.shape[1]))
        #     morphs[i] = statmorph.SourceMorphology(img_noise,seg_img,jades_id,weightmap=np.sqrt(1/wht),skybox_size=32,verbose=True)
        #     # err_ind+=1

        # #### getting error bars for concentration # just testing with y here
        # # root_asym = np.median(asym_err)
        # root_asym = root_morph.asymmetry
        # # print(root_asym)
        # asym_err = [i.asymmetry for i in morphs]
        # root_asym = np.median(asym_err)
        # lower_err = np.abs(root_asym-np.percentile(asym_err,16))
        # upper_err = np.abs(root_asym-np.percentile(asym_err,84))
        # print(f'{root_asym=}+{upper_err} & -{lower_err}')
        # plt.plot(50,root_asym,'.')
        # plt.errorbar(50, root_asym, yerr=np.array([[lower_err ,upper_err]]).T,fmt='.')
        # plt.title('asymmetry')
        # plt.ylim([0,0.8])
        # plt.xlim([0,100])
        # plt.show()

        # #### trying to get values from morphs alone
        # concs = [i.concentration for i in morphs]
        # root_conc = root_morph.concentration
        # root_conc = np.median(concs)
        # lower_conc_err = np.abs(root_conc-np.percentile(concs,16))
        # upper_conc_err = np.abs(root_conc-np.percentile(concs,84))
        # print(f'{root_conc=}+{upper_conc_err} & -{lower_conc_err}')
        # plt.plot(50,root_conc,'.')
        # plt.errorbar(50, root_conc, yerr=np.array([[lower_conc_err ,upper_conc_err]]).T,fmt='.')
        # plt.title('concentration')
        # plt.ylim([2,4.5])
        # plt.xlim([0,100])
        # plt.show()


        #### with gini & m20?
        # gini_err = [i.gini for i in morphs]
        # root_gini = root_morph.gini
        # root_gini = np.median(gini_err)
        # lower_gini_err = np.abs(root_gini-np.percentile(gini_err,16))
        # upper_gini_err = np.abs(root_gini-np.percentile(gini_err,84))
        # print(f'{root_gini=}+{upper_gini_err} & -{lower_gini_err}')
        # plt.plot(50,root_gini,'.')
        # plt.errorbar(50, root_gini, yerr=np.array([[lower_gini_err ,upper_gini_err]]).T,fmt='.')
        # plt.title('gini')
        # plt.ylim([.3,.8])
        # plt.xlim([0,100])
        # plt.show()
        
        # exit()


        ###############################################################################


        if(error_bars):
            morphs = np.zeros(100,dtype=object)
            ### try to resize the exp cutout
            # Resize the array using zoom() 
            new_size = wht.shape
            new_exp = zoom(exp, (new_size[0] / exp.shape[0], new_size[1] / exp.shape[1]), order=0) 
            exp=new_exp
            
            # plt.imshow(exp)
            # plt.show()
            # wht_with_poisson = np.where(wht.data==np.nan,0.,np.sqrt(1/wht+abs(img)/exp))

            ### following the tutorial from Gabe on how to add Poisson Noise to the WHT maps (https://dawn-cph.github.io/dja/blog/2023/07/18/image-data-products/)
            phot_scale = 1.
            phot_scale /= PHOTMJSR
            phot_scale /= PHOTSCAL
            phot_scale *= (PHOTFNU / OPHOTFNU)
            
            # print(f'source ({search_id})')
            # print(f'{PHOTMJSR=}')
            # print(f'{PHOTSCAL=}')
            # print(f'{PHOTFNU=}')
            # print(f'{OPHOTFNU=}')
            
            # effective gain e-/DN including the exposure time
            eff_gain = phot_scale*exp
            var_poisson_dn = np.maximum(img,0)/eff_gain ## max because gets rid of np.nan & negatives I think
            # this is the original vairance from the wht image = RNOISE+BG
            var_wht = 1/wht     
            var_total = var_wht+var_poisson_dn
            full_wht = 1/var_total

            # replace null weights with 0
            full_wht[var_total<=0] = 0

            # save it to fits file here if we want, otherwise just plot it
            # plt.imshow(np.sqrt(1/wht),vmin=0,vmax=0.1)
            # plt.colorbar()
            # plt.show()
            # plt.imshow(np.sqrt(1/full_wht),vmin=0,vmax=0.05)
            # plt.colorbar()
            # plt.show()
            # exit()
            # plt.imshow(wht_with_poisson,vmin=0, vmax=0.1)
            # plt.colorbar()
            # plt.show()
            # exit()

            if(use_grizli):
                # morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=wht_with_poisson,skybox_size=32,verbose=True)
                morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=np.sqrt(1/full_wht),skybox_size=32,verbose=True)
                # morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=np.sqrt(1/wht),skybox_size=32,verbose=True)
            else:
                morph = statmorph.SourceMorphology(img,seg_img,jades_id,gain=source_gain,skybox_size=32,verbose=True) # maybe cutout_extent=1?


            # fig=make_figure(morph)
            # plt.show()
            # exit()
            ### trying with the error map instead of just the snp
            # err_map = np.sqrt(1/wht)


            # print(f'{morph.concentration=}')
            # print(f'{morph.asymmetry=}')
            # print(f'{morph.gini=}')
            # print(f'{morph.m20=}')
            # exit()

            # snp = morph.sn_per_pixel
            # ### doing the scaling of the Guassian this way bc statmorph tutorial says that's how to apply bg noise to image: https://statmorph.readthedocs.io/en/latest/notebooks/tutorial.html
            # sky_sigma = 1.0 / snp # scaling of the gaussian distribution, it should ideally be based on the noise sigma of the image I think?
            fig_num=0
            for i in range(100): # size=(img.shape[0],img.shape[1])

                # img_noise = img+np.random.normal(loc=0,scale=sky_sigma,size=(img.shape[0],img.shape[1]))
                
                
                ### adding error to the image
                ## usually we want to scale the gaussian based on 1/sn
                # err = wht ## it's actually the error map, but I called it weight to make it make sense when passing into statmorph
                # err = np.sqrt(1/wht_with_poisson)
                # err = 1/wht_with_poisson**2
                # err = wht_with_poisson
                # plt.imshow(err)
                # plt.show()
                # plt.imshow(wht)
                # plt.show()
                # err = np.where(wht==0,0.,np.sqrt(1./wht))
                # err = np.sqrt(1/wht)
                # err = np.sqrt(1/wht)
                err = np.sqrt(1/full_wht) # this is the error with the poisson term
                # testing_scale = 10
                # err = 1/5 # this applies visual noise quite well but how do we convert the weightmap (values ~5000) to 0.2
                # print(np.mean(wht))
                randarr=np.random.standard_normal(img.shape)
                img_noise=img+(randarr*err)

                # snp = 100.0
                # sky_sigma = 1/snp
                # sky_sigma = np.sqrt(1/wht)
                # img_noise_2 = img+ (sky_sigma * np.random.standard_normal(size=img.shape))

                # plt.imshow(img, cmap='gray', origin='lower',norm=simple_norm(img, stretch='log', log_a=10000))
                # # plt.imshow(img, cmap='gray', origin='lower')
                # plt.show()
                # plt.imshow(img_noise, cmap='gray', origin='lower',norm=simple_norm(img_noise, stretch='log', log_a=10000))
                # plt.imshow(img_noise, cmap='gray_r', origin='lower')
                # plt.show()
                # exit()


                # # print(np.mean(img_noise_2))
                # print(np.mean(img_noise))
                # print(np.mean(img))
                # exit()

                # # err=np.where(wht.data==0,0.,np.sqrt(1./wht.data))
                # # randarr=np.random.standard_normal(sci.shape)
                # # sci+=(randarr*err)

                # print(np.mean(img))
                # print(np.mean(img_noise))
                # plt.imshow(img)
                # plt.show()

                # img_noise = img+np.random.normal(loc=0,scale=err_map,size=(img.shape[0],img.shape[1]))
                # img_noise = img+np.random.randn()*err_map*shape_somehow

                
                if(use_grizli):
                    morphs[i] = statmorph.SourceMorphology(img_noise,seg_img,jades_id,weightmap=np.sqrt(1/full_wht),skybox_size=32,verbose=False)
                    while(morphs[i].concentration==-99.0 or morphs[i].asymmetry==-99.0 or morphs[i].gini==-99.0 or morphs[i].m20==-99.0):
                    # while(morphs[i].flag>=2):
                        print(f'[!!!] Perturbation #{i} gave -99.0 for something. Rerunning with new noise...')
                        ## rerun it with new noise until it gets non error values
                        randarr=np.random.standard_normal(img.shape)
                        img_noise=img+(randarr*err)
                        morphs[i] = statmorph.SourceMorphology(img_noise,seg_img,jades_id,weightmap=np.sqrt(1/full_wht),skybox_size=32,verbose=False)
                    if(i%20==0):
                        fig = make_figure(morphs[i])
                        savepath = (f'research/statmorph_output/grizli/noise_figures/{search_id}_noise_figure_{fig_num}')
                        plt.savefig(savepath,dpi=100)
                        plt.close()
                        fig_num+=1


                else:
                    morphs[i] = statmorph.SourceMorphology(img_noise,seg_img,jades_id,gain=source_gain,skybox_size=32,verbose=True)
            
            # morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=np.sqrt(1/full_exp_cutout),skybox_size=32,verbose=True) # weight map some auto generated by SExtractor will give it as 1/RMS
            # morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=np.sqrt(1/wht),skybox_size=32,verbose=True) # weight map some auto generated by SExtractor will give it as 1/RMS
            # morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=np.sqrt(1/wht),skybox_size=32,verbose=True) # weight map some auto generated by SExtractor will give it as 1/RMS
        else:
            if(use_grizli):
                ## trying Sam's solution to add Poisson noise to the weightmap to see how it differs from without Poisson noise
                ## Poisson noise is removed only from the weightmap, not from the science image directly
                
                # print(gs_expmaps[0].shape)
                # print(gs_whtmaps[0].shape)
                # print(gs_expmaps[1].shape)
                # print(gs_whtmaps[1].shape)
                # print(gs_expmaps[2].shape)
                # print(gs_whtmaps[2].shape)
                # print(gs_expmaps[3].shape)
                # print(gs_whtmaps[3].shape)
                # print(gs_expmaps[4].shape)
                # print(gs_whtmaps[4].shape)

                ### try to resize the exp cutout
                # from scipy.ndimage import zoom 
                
                # # # Resize the array using zoom() 
                # new_size = wht.shape
                # new_exp = zoom(exp, (new_size[0] / exp.shape[0], new_size[1] / exp.shape[1]), order=0) 
                # exp=new_exp
                
                # plt.imshow(img)
                # plt.show()
                # plt.imshow(wht,origin='lower')
                # plt.show()
                # plt.imshow(exp,origin='lower')
                # plt.show()
                # plt.imshow(new_exp,origin='lower')
                # plt.show()
                
                ### instead of the line below, run the perturbations without poisson noise, just the wht map
                # wht_with_poisson = np.where(wht.data==np.nan,0.,np.sqrt(1/wht+abs(img)/exp.data))
                
                ### Getting values with Poisson noise in the weightmap
                # morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=wht_with_poisson,skybox_size=32,verbose=True)

                ### getting initial values for comparison

                morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=np.sqrt(1/wht),skybox_size=32,verbose=True)
            else:
                print('outdated method, check your input!!')
                exit()
                # morph = statmorph.SourceMorphology(img,seg_img,jades_id,gain=source_gain,skybox_size=32,verbose=True) # maybe cutout_extent=1?
        
        print(f'[---] Searching for Kirk ID: {search_id}')
        print(f'[---] Searching for JADES ID: {jades_id}')
        print(f'[---] Found JADES ID: {source["ID"][0]}')
        print(f'[>>>] Cutout size being inputted: {img.shape}')
        print(f'Source location: {position}')

        ### Get the error bar values (16th & 84th) for C, A, Gini, M20
        ### if no error bars we want the errors to be 0 for all ^^
        c_err_16,c_err_84,a_err_16,a_err_84,gini_err_16,gini_err_84,m20_err_16,m20_err_84 = [0 for i in range(8)]

        if(error_bars):
            ## morphs[i] should be set so we can get individual values

            ### concentration
            c_err = [i.concentration for i in morphs]
            # c_err = np.array([1,2,3,4,5,6,7,8,9])

            ### ensuring they're gaussian distributed
            ### so we just save all the concentration values too
            # print(f'Concentration for first source = {c_err}')
            # root_conc = np.median(c_err)
            root_conc = morph.concentration
            c_err_16 = np.abs(root_conc-np.percentile(c_err,16))
            c_err_84 = np.abs(root_conc-np.percentile(c_err,84))

            ### asymmetry
            a_err = [i.asymmetry for i in morphs]
            # a_err = np.array([1,2,3,4,5,6,7,8,9])
            # print(f'Asymmetry for first source = {a_err}')
            # root_asym = np.median(a_err)
            root_asym = morph.asymmetry
            a_err_16 = np.abs(root_asym-np.percentile(a_err,16))
            a_err_84 = np.abs(root_asym-np.percentile(a_err,84))

            ### gini
            g_err = [i.gini for i in morphs]
            # g_err = np.array([1,2,3,4,5,6,7,8,9])
            # print(f'Gini for first source = {g_err}')
            # root_gini = np.median(g_err)
            root_gini = morph.gini
            gini_err_16 = np.abs(root_gini-np.percentile(g_err,16))
            gini_err_84 = np.abs(root_gini-np.percentile(g_err,84))

            ### m20
            m_err = [i.m20 for i in morphs]
            # m_err = np.array([1,2,3,4,5,6,7,8,9])
            # print(f'M20 for first source = {m_err}')
            # root_m20 = np.median(m_err)
            root_m20 = morph.m20
            m20_err_16 = np.abs(root_m20-np.percentile(m_err,16))
            m20_err_84 = np.abs(root_m20-np.percentile(m_err,84))
        else:
            print(f'Source:{search_id}')
            print(f'Filter: {source_filter}')
            print(f'Asymmetry: {morph.asymmetry}')
            print(f'Gini:{morph.gini}')



        #### change these to be the root values of the entire distribution, not just the un-noised morph
        ########################
        row.append(source_filter)
        row.append(morph.flag)
        # row.append(morph.concentration)
        row.append(root_conc)
        row.append(c_err_16)
        row.append(c_err_84)
        row.append(str(c_err))
        # row.append(morph.asymmetry)
        row.append(root_asym)
        row.append(a_err_16)
        row.append(a_err_84)
        row.append(str(a_err))
        row.append(morph.outer_asymmetry)
        row.append(morph.smoothness)
        # row.append(morph.gini)
        row.append(root_gini)
        row.append(gini_err_16)
        row.append(gini_err_84)
        row.append(str(g_err))
        # row.append(morph.m20)
        row.append(root_m20)
        row.append(m20_err_16)
        row.append(m20_err_84)
        row.append(str(m_err))
        row.append(morph.gini_m20_merger)
        row.append(morph.multimode)
        row.append(morph.intensity)
        row.append(morph.deviation)
        row.append(morph.sn_per_pixel)
        row.append(str(img.shape))

        full_output.loc[full_output['ID']== search_id, cols] = row
        if(source_in_filter==""):
            if(id_num%round(len(search_ids)/3)==0):
                ## save the sheet every 1/3 so that if it crashes I have some output
                full_output.to_csv(f'research/statmorph_output/grizli/{output_name}.tsv','\t',index=False)
                print('Created backup file at 1/3 restore point...')




        print(f'Flag: {morph.flag}')
        # print(f'concentration \t(C): {morph.concentration} ')
        # print(f'concentration err\t {(c_err_16,c_err_84)}')
        # print(f'asymmetry (under-estimate) (A): {morph.asymmetry} ')
        # print(f'asymmetry err\t {(a_err_16,a_err_84)}')
        # print(f'outer asymmetry\t (A0): {morph.outer_asymmetry}')
        # print(f'smoothness(depends on psf heavily) (S): {morph.smoothness} ')
        # print(f'Gini: {morph.gini}')
        # print(f'gini err\t {(gini_err_16,gini_err_84)}')
        # print(f'M20: {morph.m20}')
        # print(f'm20 err\t {(m20_err_16,m20_err_84)}')
        # print(f'Gini-M20 merger stat: {morph.gini_m20_merger}')
        # print(f'Multimode (M): {morph.multimode}')
        # print(f'Intensity: {morph.intensity}')
        # print(f'Deviation: {morph.deviation}')
        print(f'S/N (per pixel): {morph.sn_per_pixel}')
        # print(f'Concentrations={str(c_err)}')
        # print(f'Asymmetries={(str(a_err))}')
        # print(f'Ginis={(str(g_err))}')
        # print(f'M20s={(str(m_err))}')

        ### for plotting the distribution of values for each, just for testing purposes
        # fig, axes = plt.subplots(2, 2, figsize=(10, 10))
        # labels = ['Concentration', 'Asymmetry', 'Gini', 'M20']
        # data = [c_err, a_err, g_err, m_err]
        # axes = axes.flatten()
        # for i in range(4):
        #     axes[i].hist(data[i], bins=20, color='b', alpha=0.7)
        #     axes[i].set_title(labels[i])
        #     axes[i].axvline(np.median(data[i]), color='r', linestyle='dashed', linewidth=1.5)
        #     fig.suptitle(f"Distribution of Values for {search_id}", fontsize=16)

        # savepath = (f'research/statmorph_output/grizli/{search_id}_Error_Histogram')
        # plt.savefig(savepath,dpi=100)
        # plt.close()
        # plt.show()
        

        # fig = make_figure(morph) # make a figure of the root morph
        savepath = (f'research/statmorph_output/{search_id}_{size}_{source_filter}')
        if(use_grizli):
            savepath = (f'research/statmorph_output/grizli/{search_id}_{size}_{source_filter}')
        # plt.savefig(savepath,dpi=100)
        # plt.show()
        # plt.close()
        # exit()
        id_num+=1

    full_output.to_csv(f'research/statmorph_output/grizli/{output_name}.tsv','\t',index=False)
    print('[!!!] Finished statmorph measurements -- outputing table...')
    print(full_output[['ID','AGN^a','z^b','Filter','Flag','Concentration (C)','Asymmetry (A)','Gini','M20','S/N','Cutout Size (pix)']])
    disclaimer()
    return

##################

def test_in_ceers(note=''):
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
    # search_one_id = np.array([11689])
    # for id in search_one_id:
    for id in search_ids:
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

        print('------------------------------')
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
        # plt.imshow(deblended_segmap,origin='lower')
        # plt.show()
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
        morph_full = statmorph.source_morphology(img,deblended_segmap,gain=source_gain,cutout_extent=1.5,verbose=True) # maybe cutout_extent=1?
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
        stat_rows[c,12] = (f'{img.shape}') # size of cutout
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
        
        # fig = make_figure(morph)
        # plt.show()
        # plt.savefig(f'/home/robbler/research/CEERS_statmorph_testing/statmorph_fits/{id}.png')
        # plt.close()
    statmorph_table = Table(rows=stat_rows, names=['id', 'position','flag','concentration (C)','asymmetry (A)','outer asymmetry (Ao)','deviation','smoothness (S)','gini','m20','gini/m20 merger','s/n ratio','cutout size'])
    table_output = statmorph_table.to_pandas()
    print(table_output)
    table_output.to_csv(f"/home/robbler/research/CEERS_statmorph_testing/{note}-CEERS_statmorph.csv", encoding='utf-8', index=False)



# print('with inputting the psf -- no adjustments base sample test')
# test_in_ceers(note='without-psf')
# exit()

# grizli_filters = ['f115w','f150w','f200w','f277w','f356w','f444w']
# grizli_filters = ['f277w','f356w','f444w'] # testing with these bc others are giving errors?
# grizli_filters = ['f356w']
# jades_filters = ['f115w','f150w','f200w','f277w','f356w','f444w']


# sources = ['GN_IRS11','GN_IRS17','GN_IRS27','GS_IRS14','GS_IRS25','GS_IRS60','GS_IRS9']
filters = ['f150w','f200w','f277w','f356w','f444w']
# for f in filters:
#     try:
run_in_jades(use_grizli=True,error_bars=True,output_name='statmorph-non277-rest-12-23',source_in_filter='f277w')
    # except:
    #     continue




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
# command 2> >(tee /home/robbler/research/statmorph_output/logs/$(date +"%Y_%m_%d_%I_%M_%p")_err) 1> >(tee /home/robbler/research/statmorph_output/logs/$(date +"%Y_%m_%d_%I_%M_%p")_out) | tee >/home/robbler/research/statmorph_output/logs/$(date +"%Y_%m_%d_%I_%M_%p")_all.txt