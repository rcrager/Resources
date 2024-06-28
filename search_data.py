from astropy.table import Table
from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u



############ README #################
'''
Script that has useful functions for finding objects, or IDs inside fits files
'''
#####################################


# general structure, need to edit this later...
def get_subcatalog(filename,search_ids, args=('id','ra','dec')):
    '''
    Function to search for a list of object IDs and parameters (usually ra,dec) from a larger catalog of data.
    '''
    catalog_filename = '/home/robbler/research/CEERS_statmorph_testing/CEERS_PR1_LW_SUPER_catalog.fits'
    
    with fits.open(catalog_filename) as hdul:
        # this way it closes automatically when with: ends
        flag_data = Table(hdul[1].data)
        table = flag_data['id','x','y','ra','dec']
        search_ids = np.array([13447,21320,19427,6944,6611,31265,19511,13206,2705,14053,20998,918,31485,3747,11689,20248,11790,23584,23419,13219,16315,12403,2985])
        c=0
        rows = np.zeros((len(search_ids),3))
        for i in search_ids:
            mask = table['id']==i
            rows[c,0] = int(i)
            rows[c,1] = table[mask]['ra']
            rows[c,2] = table[mask]['dec']
            c+=1
        table = Table(rows=rows, names=['id', 'ra','dec'])
        return table
    


########### edit this later for more general use?  ##############
###########, more so keeping it here for reference ##############
def get_obj_near(filename,search_ids,size=3,args=('id','ra','dec')):
    '''
    Function to pull nearest neighbors within a certain radius (size in arcseconds) of a given set of ids (looking in catalog given by filename)
    '''
    # opens fits file (default to read only)
    fits_image_filename = 'research/JADES catalog/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits'
    with fits.open(fits_image_filename) as hdul:
        tolerance = size
        hdr = hdul[0].header
        flag_data = Table(hdul[2].data)
        
        catalog = SkyCoord(flag_data['RA']*u.deg,flag_data['DEC']*u.deg)
        
        # kirk_galaxy_coords comes from the search_ids in the catalog, so we'd need to loop through the catalog to find them 
        search = SkyCoord(kirk_galaxy_coords[:,0]*u.deg,kirk_galaxy_coords[:,1]*u.deg)

        ### match two catalogs of data
        ### gets all matches for the search catalog to the full catalog within the tolerance
        ### then we can just loop through this catalog[idxc] list to create the regions file
        ### https://keflavich-astropy.readthedocs.io/en/latest/coordinates/matchsep.html
        idx_search, idx_catalog, d2d, d3d = catalog.search_around_sky(search, tolerance*u.arcsec)  
        print(len(catalog[idx_catalog]))
        coords_match = catalog[idx_catalog]


def make_regions(table,radius=1):
    ##### doesn't look like this is a good idea, for now it's so much easier to print to output then save to file using manual region texts
    ##### then we can just import from regions file to get fits file or so from here: https://astropy-regions.readthedocs.io/en/stable/region_io.html
    from regions import Region, CircleSkyRegion
    '''
    function to make ds9 regions file (circles) from an astropy table, table should ideally be in ASCII format so it's easier for me
    this is just meant to see the subcatalog you've created; radius in arcseconds
    '''
    for i in table:
        coord = SkyCoord(i['ra'], i['dec'], unit='deg', frame='fk5')
        region = CircleSkyRegion(coord, radius=radius * u.arcsec)
        region.meta['label'] = i['id']
        print(region.serialize(format='ds9'))
     #   region.write('/home/robbler/research/Scripts/test_regions_file.reg',format='ds9',overwrite=True)
    # coord = SkyCoord(table['ra'],table['dec'],unit='deg',frame='fk5')
    # regions = CircleSkyRegion(coord,radius=radius*u.arcsec)
    # regions.write('/home/robbler/research/Scripts/test_regions_file.reg',format='ds9',overwrite=True)
    #print(coord)