
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib
from astropy.visualization import simple_norm
from astropy.modeling.models import Sersic2D
from astropy.convolution import convolve, Gaussian2DKernel, CustomKernel
import pandas as pd
from scipy.ndimage import zoom

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
import statmorph
from photutils.segmentation import make_2dgaussian_kernel, SegmentationImage
from astropy.stats import SigmaClip
#sigma_clipping = SigmaClip(sigma=sigma)

# for finding the psf
# import pypher

from astropy.utils.exceptions import AstropyWarning


def do_morph():
    ######## store everything in respective array location for the following Nircam filters
    nir_filters = ['f150w','f200w','f277w','f356w','f444w']

    psf_locs = [f'/home/robbler/research/JADES catalog/PSF_NIRCam_in_flight_opd_filter_{i.upper()}.fits' for i in nir_filters]

    psfs = [fits.open(f) for f in psf_locs]
    psf_data = [f[0].data for f in psfs]
    psf_pixscl =  [f[0].header['PIXELSCL'] for f in psfs]
    [f.close() for f in psfs]

    conv_nir_filters = ['f150w','f200w','f277w','f356w']

    ### now open the combined psf kernels to later convolve them with the images to get them to f277
    kernel_psf_locs = [f'/home/robbler/research/JADES catalog/{i}_to_f444w_fit_pix_scale.fits' for i in conv_nir_filters]
    psf_kernel_to_444 = [fits.getdata(f) for f in kernel_psf_locs]


    #### need to open f356 & f444 to get the following sources measurement differences for C,A,Gini,M20
    # GN_IRS36
    # GN_IRS43
    # GS_IRS9
    test_ids = ['GN_IRS36','GN_IRS43','GS_IRS9','GN_IRS4']



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

    # maybe this is the weight value we're looking for for statmorph (one value for each entry, so how would I convert this to an image to be used with statmorph though?)
    #print(jades_catalog[3]['F444W_WHT'])

    # fwhms = {'F070W':0.742,'F090W':0.968,'F115W':1.194,'F150W':1.581,'F200W':2.065,'F277W':1.397,'F356W':1.810,'F444W':2.222}
    # # ^ as described in https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-point-spread-functions#gsc.tab=0

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
    # search_ids = test_ids

    search_ids = ['GS_IRS12',
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
    'GS_IRS9'] 
    

    ### search ids for testing error bar distribution
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




    cols = ['Flags','C_f150w','C_f200w','C_f277w','C_f356w','C_f444w','A_f150w','A_f200w','A_f277w','A_f356w','A_f444w','Gini_f150w','Gini_f200w','Gini_f277w','Gini_f356w','Gini_f444w','M20_f150w','M20_f200w','M20_f277w','M20_f356w','M20_f444w','ARMS_f150w','ARMS_f200w','ARMS_f277w','ARMS_f356w','ARMS_f444w']
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

        data_source = sources[(sources['ID'].str.contains(id,regex=False))].iloc[0]
        source_row = sources.loc[sources['ID'] == id]
        source_filter = source_row['Obs Filter'].iloc[0]
        
        flags = np.zeros(len(nir_filters))
        cfilt = np.zeros(len(nir_filters))
        afilt = np.zeros(len(nir_filters))
        ginifilt = np.zeros(len(nir_filters))
        m20filt = np.zeros(len(nir_filters))
        armsfilt = np.zeros(len(nir_filters))
        for source_filter in nir_filters: # loop through all filters and get the c,a,gini,m20 then write to tsv
            # tmp_filters = ['f150w','f200w','f277w','f356w','f444w']
            cur_psf = psf_data[nir_filters.index(source_filter)]


            # index = conv_nir_filters.index(source_filter)
            index = nir_filters.index(source_filter)

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
                    # seg = Cutout2D(jades_seg,position,0.75*size*units.arcsec,wcs=jades_wcs_coords,copy=True).data

                # #### testing full_exp for weight map instead because it includes the source poisson noise?
                # full_exp_cutout = Cutout2D(full_exp_img,position,size*units.arcsec,wcs=wcs_coords,copy=True).data
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

            # if(source_filter!='f444w'):
            #     tmp_kern = psf_kernel_to_444[nir_filters.index(source_filter)]
            #     target_kern = tmp_kern[0:-1,0:-1] # to make it odd axes because for some reason astropy will NOT LET YOU DO EVEN AXES AHHHH WHY


            #     # print(f'{tmp_kern.shape=}')
            #     # print(f'{target_kern.shape=}')

            #     # fig,ax = plt.subplots(2,1)
            #     # ax[0].imshow(tmp_kern)
            #     # ax[1].imshow(target_kern)
            #     # plt.show()


            #     psf_kernel = CustomKernel(target_kern)
            #     # print(img.shape)
            #     # print(psf_kernel.shape)
            #     conv = convolve(img,psf_kernel)
            #     ### astropy.convolution.utils.KernelSizeError: Kernel size must be odd in all axes.
            #     ### now run statmorph on this with error bars & see what happens

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

            def poisson(image):
                var_poisson_dn = np.maximum(image,0)/eff_gain ## max because gets rid of np.nan & negatives I think
                # this is the original vairance from the wht image = RNOISE+BG
                var_wht = 1/wht     
                var_total = var_wht+var_poisson_dn
                full_wht = 1/var_total

                # replace null weights with 0
                full_wht[var_total<=0] = 0
                return full_wht
            
            full_wht = poisson(img)
            # if(source_filter=='f444w'):
            #     full_wht_conv = full_wht
            #     conv = img
            # else:
            #     full_wht_conv = poisson(conv)

            # morph = statmorph.SourceMorphology(conv,seg_img,jades_id,weightmap=np.sqrt(1/full_wht),psf=cur_psf,skybox_size=32,verbose=True)
            # morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=np.sqrt(1/full_wht),psf=cur_psf,skybox_size=32,verbose=True)
            # morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=np.sqrt(1/full_wht),skybox_size=32,verbose=True)
            morph = statmorph.SourceMorphology(img,seg_img,jades_id,weightmap=np.sqrt(1/full_wht),skybox_size=32,verbose=False)
            # conv_morph = statmorph.SourceMorphology(conv,seg_img,jades_id,weightmap=np.sqrt(1/full_wht_conv),skybox_size=32,verbose=True)

            
            # row.append(source_filter)
            # row.append(morph.flag)
            # row.append(conv_morph.concentration)
            flags[nir_filters.index(source_filter)] = morph.flag
            cfilt[nir_filters.index(source_filter)] = morph.concentration
            afilt[nir_filters.index(source_filter)] = morph.asymmetry
            ginifilt[nir_filters.index(source_filter)] = morph.gini
            m20filt[nir_filters.index(source_filter)] = morph.m20
            armsfilt[nir_filters.index(source_filter)] = morph.rms_asymmetry2

            fig = make_figure(morph)
            plt.savefig(f'research/making_psf_edits/k-correction-make-figures/new_segmap_{search_id}_{source_filter}.png',dpi=100)
            plt.close()
            # row.append(morph.concentration)
            # # row.append(conv_morph.asymmetry)
            # row.append(morph.asymmetry)
            # # row.append(conv_morph.gini)
            # row.append(morph.gini)
            # # row.append(conv_morph.m20)
            # row.append(morph.m20)
            # # row.append(conv_morph.sn_per_pixel)
            # row.append(str(img.shape))
        


        print(f'Done with {search_id}...')
        row.append(str(flags))
        [row.append(c) for c in cfilt]
        [row.append(c) for c in afilt]
        [row.append(c) for c in ginifilt]
        [row.append(c) for c in m20filt]
        [row.append(c) for c in armsfilt]


        full_output.loc[full_output['ID']== search_id, cols] = row
        # if(id_num%round(len(search_ids)/4)==0):
            ## save the sheet every 1/4 so that if it crashes I have some output
            # full_output.to_csv(f'research/making_psf_edits/k-correction-measurements.tsv','\t',index=False)
            # print('Created backup file at 1/4 restore point...')
        id_num+=1


    full_output.to_csv(f'research/making_psf_edits/gs_k-correction-measurements.tsv','\t',index=False)
    print('[!!!] Finished statmorph measurements -- outputing table...')
    print(full_output[['ID','AGN^a','Flags','ARMS_f150w','ARMS_f200w','ARMS_f277w','ARMS_f356w','ARMS_f444w']])
    disclaimer()


################################################
### for graphing purposes
relpath = 'research/'
small_size = 16
large_size = 18
matplotlib.rc('font', size=large_size)
matplotlib.rc('axes', titlesize=large_size)
matplotlib.rc('legend',fontsize=small_size)
#################################################

def plot_conv(output,search,name): # since then I've edited the script to look for k-correction instead...
    ### make plots comparing the convolved with the original
    agn = output['AGN^a']
    og = output[f'{search} Original']
    conv = output[f'{search} Convolved']
    sf = np.where(agn<20)[0]
    comp = np.where((agn>=20)&(agn<80))[0]
    agns = np.where(agn>=80)[0]

    x_key = ['Original','Convolved']
    cols = ['green','orange','blue']
    names = ['SF','Composite','AGN']

    perc_diff = np.abs(og-conv)/og
    plt.scatter(agn,perc_diff,marker='o')
    # c=0
    # for ind, row in output.iterrows():
    #     if(row['C_f150w']!=np.nan):
    #         if(row['AGN^a']<20):
    #             c=0
    #         elif((row['AGN^a']>=20)&(row['AGN^a']<80)):
    #             c=1
    #         else: #it's agn >=80
    #             c=2
                
    #             y_vals = np.array([np.abs((obs_val-row[f'{search}_{f}'])/obs_val) for f in filters])
    #             y_vals*=100
    #             # zero = y_vals[y_vals == 0]
    #             x_ind = filters.index(obs_filt)
    #             line = plt.plot(x_key,y_vals,color=cols[c],marker='o',zorder = 0)
    #             # print(line[0]._y[line[0]._y == 0.])
    #             scat = plt.scatter(x_key[x_ind],y_vals[x_ind],marker='o',color='red',zorder=1)
    #             line[0].set_label(names[c])
    # for ind in [sf,comp,agns]:
    #     xs = ([0 for i in range(len(og[ind]))],[1 for i in range(len(og[ind]))])
    #     line = plt.plot(xs,(og[ind],conv[ind]),color=cols[c],marker='o')
    #     line[0].set_label(names[c])
    #     c+=1

    # sc = plt.scatter(xs[0],conc_og,c=agn,cmap='viridis',linestyle='-',marker='o')
    # plt.scatter(xs[1],conc_conv,c=agn,cmap='viridis',linestyle='-',marker='^')
    # fig.canvas.draw()
    # face_colors = sc.get_facecolors()
    # col = sc.get_facecolors()[-2]
    # ## now draw the lines
    # print(col)
    # plt.plot(x_key,(conc_og,conc_conv),color=col)

    plt.ylabel(f'{name} (% Diff)')
    plt.xlabel('AGN(%)')
    # plt.xticks(ticks=[0, 1], labels=x_key)
    plt.legend()

    # sc = ax.scatter(x, y, c=y, cmap=cmap, norm=norm) # Add the color bar 
    # colorbar = plt.colorbar() 
    # colorbar.set_label('$f_{{AGN}}$')
    # plt.colorbar()
    plt.savefig(f'{relpath}making_psf_edits/{name}_plot.png',dpi=200)
    plt.show()

def plot_k_corr(output,search,ax_name,perc_diff=False,plot_diff=False):
    # plot the % diff [(original-new)/original] values for each filter (for each source) on one plot here (I think)
    # maybe sort by agn fraction again and color code into the 3 categories
    # agn = output['AGN^a']
    filters = ['f150w','f200w','f277w','f356w','f444w']
    # y_s = output[[f'{search}_f150w',f'{search}_f200w',f'{search}_f277w',f'{search}_f356w',f'{search}_f444w']]
    # obs = output['Obs Filter']

    x_key = [1.501,1.988,2.776,3.566,4.401]

    cols = ['green','orange','blue']
    names = ['SF','Composite','AGN']

    scatter_vals = np.full(len(output),fill_value=0.0)
    s=0
    for ind, row in output.iterrows():
        # if(row['C_f150w']!=0):
        if(row['AGN^a']<20):
            c=0
        elif((row['AGN^a']>=20)&(row['AGN^a']<80)):
            c=1
        else: #it's agn >=80
            c=2
            # print(f'AGN name: {row["ID"]} with asym {row["A_f277w"]}')
        obs_filt = row['Obs Filter']
        obs_val = row[f'{search}_{obs_filt}']
        # y_vals = np.array([np.abs((obs_val-row[f'{search}_{f}'])/obs_val) for f in filters])
        # y_vals*=100
        # zero = y_vals[y_vals == 0]
        x_ind = filters.index(obs_filt)
        targ_filt = 'f277w'
        targ_ind = filters.index(targ_filt)
        if(plot_diff):
            ## plot the new asymmetry vs old asymmetry
            targ_val = row[f'{search}_{targ_filt}']
            scatter_vals[s] = targ_val-obs_val # not doing abs here because it's higher error this way (not sure why but statistics stuff probably)
            line = plt.errorbar(targ_val,obs_val,yerr=scatter_vals[s],color=cols[c],marker='o',zorder = 0,label=names[c])
            # err = plt.errorbar(targ_val,obs_val,yerr=scatter_vals[c],marker='o',color='red')
        else:
            if(targ_ind<x_ind): 
                new_filts = filters[targ_ind:x_ind+1]
                new_x_key = x_key[targ_ind:x_ind+1]
            elif(targ_ind>x_ind): 
                new_filts = filters[x_ind:targ_ind+1]
                new_x_key = x_key[x_ind:targ_ind+1]
            else: 
                new_filts = filters[x_ind:x_ind+1] # if they're = just one dot
                new_x_key = x_key[x_ind:x_ind+1]

        # print(new_filts)
        # if(perc_diff):
        #     # y_vals = np.array([np.abs(obs_val-row[f'{search}_{f}']/obs_val) for f in new_filts])
        #     y_vals = np.array([np.abs(obs_val-row[f'{search}_{f}'])/obs_val for f in new_filts])
        # else:
        #     continue
            # y_vals = np.array([np.abs((row[f'{search}_{f}'])/obs_val) for f in new_filts])
        # if((y_vals>1.1).any() or (y_vals<0.9).any()):
        #     print(f'Source Outside Tolerance: {row["ID"]}-{search}')
        # if(len(new_x_key)!=1): # if it's already in the 277 then just plot red dot there
        # line = plt.plot(new_x_key,y_vals,color=cols[c],marker='o',zorder = 0)
        # line = plt.plot(new_x_key,y_vals,color=cols[c],marker='o',zorder = 0)

        # line[0].set_label(names[c])
        # print(line[0]._y[line[0]._y == 0.])
        # x_filt_ind = new_filts.index(obs_filt)
        # scat = plt.scatter(new_x_key[x_filt_ind],y_vals[x_filt_ind],marker='o',color='red',zorder=1)
        s+=1
    # plt.ylabel(f'{ax_name} (New/Original)')
    

    plt.ylabel(f'$\sim 1 \mu m$ {ax_name}')
    plt.xlabel(f'f277w {ax_name}')
    plt.grid(True)
    x = np.linspace(np.min(output[f'{search}_f277w']), np.max(output[f'{search}_f277w']), 50) 
    y = x 
    plt.plot(x,y,color='black',linestyle='--')

    ## find scatter & plot error bars
    for s in scatter_vals:
        print(f'{s},',end='')
    std = np.std(scatter_vals)
    print(f'for {ax_name} standard error = {std}')
    plt.plot(x,y+std,color='red',linestyle=':')
    plt.plot(x,y-std,color='red',linestyle=':')

    # plt.xlim(1.9,3.6)
    # plt.ylim(0,30)
    # plt.axhline(1.1,linestyle=':',color='gray',zorder=2)
    # plt.axhline(0.9,linestyle=':',color='gray',zorder=2)
    # plt.axhline(0.1,linestyle=':',color='gray',zorder=2)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    # plt.legend()
    plt.tight_layout(pad=0.5)
    # plt.savefig(f'{relpath}making_psf_edits/newnew__new_{search}_vs_old_{search}.png',dpi=200)
    # plt.savefig(f'{relpath}making_psf_edits/_k_correction_{ax_name}.png',dpi=200)

    plt.show()

    # obs_val = output[f'{search}_{obs}']
    # # y_s = output[[f['{search}_{f}'] for f in filters]]
    # # y_s = []
    # # for i in filters:
    # #     y_s.append(output[f'{search}_{i}'])
    # # og = output[f'{search}_{filt}']
    # sf = np.where(agn<20)[0]
    # comp = np.where((agn>=20)&(agn<80))[0]
    # agns = np.where(agn>=80)[0]


    # c=0
    # for ind in [sf,comp,agns]:
    #     y_vals = [(obs_val-y_s[f'{search}_{f}'][ind])/obs_val for f in filters]
    #     line = plt.plot((x_key),y_vals,color=cols[c],marker='o')
    #     line[0].set_label(names[c])

    #     c+=1
    # plt.ylabel(ax_name)
    # plt.legend()
    # plt.savefig(f'{relpath}making_psf_edits/k_correction_{ax_name}.png',dpi=200)
    # plt.show()

    # return

#### need to edit to make the old k_correction_asymmetry figure
def plot_old_k(output,search,ax_name):
    # plot the % diff [(original-new)/original] values for each filter (for each source) on one plot here (I think)
    # maybe sort by agn fraction again and color code into the 3 categories
    # agn = output['AGN^a']
    filters = ['f150w','f200w','f277w','f356w','f444w']
    # y_s = output[[f'{search}_f150w',f'{search}_f200w',f'{search}_f277w',f'{search}_f356w',f'{search}_f444w']]
    # obs = output['Obs Filter']

    x_key = [1.501,1.988,2.776,3.566,4.401]

    cols = ['green','orange','blue']
    names = ['SF','Composite','AGN']

    for ind, row in output.iterrows():
        # if(row['C_f150w']!=0):
        if(row['AGN^a']<20):
            c=0
        elif((row['AGN^a']>=20)&(row['AGN^a']<80)):
            c=1
        else: #it's agn >=80
            c=2
        obs_filt = row['Obs Filter']
        obs_val = row[f'{search}_{obs_filt}']
        # y_vals = np.array([np.abs((obs_val-row[f'{search}_{f}'])/obs_val) for f in filters])
        # y_vals*=100
        # zero = y_vals[y_vals == 0]
        x_ind = filters.index(obs_filt)
        targ_filt = 'f277w'
        targ_ind = filters.index(targ_filt)
        if(targ_ind<x_ind): 
            new_filts = filters[targ_ind:x_ind+1]
            new_x_key = x_key[targ_ind:x_ind+1]
        elif(targ_ind>x_ind): 
            new_filts = filters[x_ind:targ_ind+1]
            new_x_key = x_key[x_ind:targ_ind+1]
        else: 
            new_filts = filters[x_ind:x_ind+1] # if they're = just one dot
            new_x_key = x_key[x_ind:x_ind+1]
        # new_filts = filters
        # new_x_key = x_key
        print(f'Source: {row["ID"]} is {names[c]} and going through filters: {new_filts}')
        y_vals = np.array([np.abs(obs_val-row[f'{search}_{f}'])/obs_val for f in new_filts])
        # y_vals*=100
        # if((y_vals>1.1).any() or (y_vals<0.9).any()):
        #     print(f'Source Outside Tolerance: {row["ID"]}-{search}')
        # if(len(new_x_key)!=1): # if it's already in the 277 then just plot red dot there
        # line = plt.plot(new_x_key,y_vals,color=cols[c],marker='o',zorder = 0)
        line = plt.plot(new_x_key,y_vals,color=cols[c],marker='o',zorder = 0)

        line[0].set_label(names[c])
        # print(line[0]._y[line[0]._y == 0.])
        x_filt_ind = new_filts.index(obs_filt)
        scat = plt.scatter(new_x_key[x_filt_ind],y_vals[x_filt_ind],marker='o',color='red',zorder=1)
    plt.ylabel(f'{ax_name} (% diff)')
    # plt.ylabel(f'Rest-Frame {ax_name}')
    # plt.xlabel(f'f277w {ax_name}')
    plt.grid(True)

    # plt.xlim(1.9,3.6)
    plt.ylim(0,.55)
    # plt.axhline(1.1,linestyle=':',color='gray',zorder=2)
    # plt.axhline(0.9,linestyle=':',color='gray',zorder=2)
    # plt.axhline(0.1,linestyle=':',color='gray',zorder=2)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    # plt.legend()
    # plt.savefig(f'{relpath}making_psf_edits/new_{search}_vs_old_{search}.png',dpi=200)
    plt.savefig(f'{relpath}making_psf_edits/_old_k_correction_{ax_name}.png',dpi=200)
    plt.show()



# do_morph()

# output = pd.read_csv("{relpath}making_psf_edits/statmorph-convolved-meas.tsv",sep='\t')
# plot_conv(output,'C','Concentration')
# plot_conv(output,'A','Asymmetry')
# plot_conv(output,'Gini','Gini')
# plot_conv(output,'M20','M20')


in_data = pd.read_csv(f"{relpath}making_psf_edits/k-correction-measurements.tsv",sep='\t')

def plot_by_agn(indata, search, ax_name): #plot agn on x-axis & measurement on y
    filters = ['f150w','f200w','f277w','f356w','f444w']
    cols = ['green','orange','blue']
    for ind,row in indata.iterrows():
        if(row['AGN^a']>=80): c=2
        if((row['AGN^a']<80)&(row['AGN^a']>=20)): c=1
        if(row['AGN^a']<20): c=0
        x_ind = filters.index(row['Obs Filter'])
        targ_filt = 'f277w'
        # if(row[f'{search}_{targ_filt}']==0):continue
        targ_ind = filters.index(targ_filt)
        if(targ_ind<x_ind): 
            new_filts = filters[targ_ind:x_ind+1]
        elif(targ_ind>x_ind): 
            new_filts = filters[x_ind:targ_ind+1]
        else: 
            new_filts = filters[x_ind:x_ind+1] # if they're = just one dot
        y_err = np.std(np.array([row[f'{search}_{fi}'] for fi in new_filts]))

        plt.errorbar(row['AGN^a'],row[f'{search}_f277w'],marker='o',yerr=y_err,color=cols[c])
    plt.ylabel(f'{ax_name}')
    # plt.ylim(-.1,1.0)
    plt.xlabel("AGN(%)")
    plt.savefig(f'{relpath}making_psf_edits/agn_{search}_plot.png')
    plt.show()
# plot_by_agn(in_data, 'A','Asymmetry')
# plot_by_agn(in_data, 'C','C')
# plot_by_agn(in_data, 'Gini','Gini')
# plot_by_agn(in_data, 'M20','M20')

# plot_old_k(in_data,'A','Asymmetry')

# plot_k_corr(in_data,'ARMS','$A_{{RMS}}$')

plot_k_corr(in_data,'A','Asymmetry',plot_diff=True)

plot_k_corr(in_data,'C','Concentration',plot_diff=True)
plot_k_corr(in_data,'Gini','Gini',plot_diff=True)
plot_k_corr(in_data,'M20','M20',plot_diff=True)