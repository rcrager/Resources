from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u
import numpy as np
from import_data import get_sources

# before opening the fits file, get the new coordinate ranges to look for
# search for each ra and dec within a certain range (+-1 deg) of the Kirkpatrick source coordinates



#### format is [RA_1, DECL_1],[ra_2,decl_2],...
#kirk_galaxy_coords = np.fromstring(raw_GN_data,sep=' ').reshape((-1,2))
#kirk_galaxy_names = raw_GN_names.strip().split('\n')

source_data = get_sources()

# tolerance to search for (arcsec) (will be +-)
tolerance = 2.0
'''
### for the GN format is now swapped to [:,0] for RA, vs [0,:] for ra from GS
kirk_low_ra = kirk_galaxy_coords[:,0]-tolerance
kirk_low_dec = kirk_galaxy_coords[:,1]-tolerance 

kirk_high_ra = kirk_galaxy_coords[:,0]+tolerance
kirk_high_dec = kirk_galaxy_coords[:,1]+tolerance
'''

### any time we make a small adjustment (+ or -) to the coordinates, need to add the cos(decl.) factor to the RA (no additional change to decl)
#kirk_low_ra = (kirk_galaxy_coords[:,0]-tolerance)*np.cos(np.radians(kirk_galaxy_coords[:,1]))
#kirk_high_ra = (kirk_galaxy_coords[:,0]+tolerance)*np.cos(np.radians(kirk_galaxy_coords[:,1]))


# make an np array 10 long for each source to allow 10 jades matches for each source 
    # (if any of them meet the max 10 source matches, then we can up the limit)
coords_match = np.zeros([10,3,10])
# array structure is matched with index of galaxy name
# so i=0 = kirk_galaxy_names[0] = kirk_galaxy_coords[0,0] = kirk_galaxy_coords[1,0]
# then for this coords_match[:,0] = JADES ID
# & coords_match[:,1] = ra
# & coords_match[:,2] = dec

# opens fits file (default to read only)
fits_image_filename = 'research/JADES catalog/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits'
with fits.open(fits_image_filename) as hdul:
    # this way it closes automatically when with: ends
    hdr = hdul[0].header
    flag_data = Table(hdul[2].data)
    size_data = hdul[3].data
    search_table = flag_data['ID','RA','DEC']
    
    catalog = SkyCoord(flag_data['RA']*u.deg,flag_data['DEC']*u.deg)
    #tolerance = 1.0 # in arcsec
    search = SkyCoord(kirk_galaxy_coords[:,0]*u.deg,kirk_galaxy_coords[:,1]*u.deg)

    ### match two catalogs of data
    ### gets all matches for the search catalog to the full catalog within the tolerance
    ### then we can just loop through this catalog[idxc] list to create the regions file
    ### https://keflavich-astropy.readthedocs.io/en/latest/coordinates/matchsep.html
    idx_search, idx_catalog, d2d, d3d = catalog.search_around_sky(search, tolerance*u.arcsec)  
    print(len(catalog[idx_catalog]))
    coords_match = catalog[idx_catalog]
    

    '''
    i = 0
    for ra_low in kirk_low_ra:
        indices = np.where((search_table['RA']>kirk_low_ra[i]) & (search_table['RA']<kirk_high_ra[i]) &
                              (search_table['DEC']>kirk_low_dec[i]) & (search_table['DEC']<kirk_high_dec[i]))[0]
        #jades_coords = np.zeros([1,3,10])
        jades_coords = np.array([search_table[indices]['ID'],search_table[indices]['RA'],search_table[indices]['DEC']])
        # example 
#        arr = np.array([data['ra'], data['dec']])
        #print(jades_coords)
        coords_match[i] = np.pad(jades_coords,((0,0),(0,len(coords_match[i][0])-len(jades_coords[0]))), mode='constant')
        #print(coords_match[i])
        if(i==21): # looking at only the GN_IRS55 HERE
            print(f'{kirk_low_ra[i]=}')
            print(f'{kirk_high_ra[i]=}')
            print(f'{kirk_low_dec[i]=}')
            print(f'{kirk_high_dec[i]=}')

            
            #print(undetected)
            # print the tolerance and bounds of ra and dec here, see if they containt the location of jades ID:1027288
            print(jades_coords)
            print(coords_match[i])


        #jades_coords.resize(3,10)
        
        #coords_match[i] = jades_coords
        #print(kirk_low_ra[i]<kirk_high_ra[i])
        #print((search_table['RA']<kirk_high_ra[i]))

        #new_data = hdul[2].data
        #m_tmp = ((new_data['RA'] > kirk_low_ra[i]) & (new_data['RA']<kirk_high_ra[i]))
        #print(m_tmp)
        
        #data = hdul[1].data
        #mask = data['mag'] > -0.5
        #newdata = data[mask]
        i+=1
        '''

### now print the new regions file
def make_regions(ra,dec,label,size,title=None):
    '''function to make regions file given equal size ra,dec,label arrays with size radius.
        ra & dec can be in (") or (deg)
        label is string
        radius should be in (")'''
    c = 0
    for i in ra:
        if ra[c]==0.0:
            break
        print(f'point({ra[c]}, {dec[c]})  # point=x {size} text={{{int(label[c])}}}')
        c+=1

def make_circles(ra,dec,label):
    '''function to make regions file given equal size ra,dec,label arrays with size radius.
        ra & dec can be in (") or (deg)
        label is string
        radius should be in (")'''
    c = 0
    for i in ra:
        if ra[c]==0.0:
            break
        print(f'circle({ra[c]}, {dec[c]}, 1.0")  # text={{{label[c]}}}')
        c+=1

pt_size = 20


# init statement

print("""# Region file format: DS9 version 4.1
global color=cyan dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5""")
    
print('# JADES Coordinates')
#print regions for jades coords
for i in coords_match:
    print(f'point({i.ra.deg}, {i.dec.deg})  # point=x {pt_size}')
    #make_regions(i[1],i[2],i[0],pt_size,title='JADES Coordinates')
    #make_regions(i.ra.deg,i.dec.deg,0,pt_size,title='JADES Coordinates')

# print regions for kirk circle coordinates
print("""global color=black""")
print('# Kirkpatrick Coordinates')
c=0
make_circles(kirk_galaxy_coords[:,0],kirk_galaxy_coords[:,1],kirk_galaxy_names)