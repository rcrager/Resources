import numpy as np
from PIL import Image, ImageOps
from import_data import get_GS,get_GN

#### Combine the images to make the form images

### get the source names & jades ids from spreadsheet
#### most recently updated for GN data (directory names n stuff)
sources = get_GN()
kirk_names = sources[0]
jades_names = sources[1]

c=0
for i in kirk_names:
    ## loop through the galaxy name values referencing the jades values too (need for jades viewer images)
    path_1 = (f'research/Galaxies/cutouts/GN/{kirk_names[c]}_f150w.png')
    path_2 = (f'research/Galaxies/cutouts/GN/{kirk_names[c]}_f277w.png')
    path_3 = (f'/home/robbler/research/Galaxies/JADES_viewer_images/GN/{jades_names[c]}.png')
    list_im = [path_1,path_2,path_3]
    imgs = [Image.open(i) for i in list_im]

    padd = 75
    imgs[0] = imgs[0].crop((0,padd,imgs[0].width,imgs[0].height-padd))
    imgs[1] = imgs[1].crop((0,padd,imgs[1].width,imgs[1].height-padd))

    '''
    dst = Image.new('RGB', (imgs[0].width+imgs[1].width +imgs[2].width, imgs[0].height),(0,255,255))

    dst.paste(imgs[0], (0, 0))
    dst.paste(imgs[1], (imgs[0].width, 0))
    dst.paste(imgs[2], (imgs[0].width+imgs[1].width, 0))
    filename = f'{kirk_names[c]}_{jades_names[c]}.jpg'
    dst.save(f'/home/robbler/research/Galaxies/form_images/{filename}')
    print(f'Saved file: {filename}')
    exit()
    '''
    #dst = Image.new('RGB', (imgs[0].width + imgs[1].width + imgs[2].width, max(imgs[0].height,imgs[1].height,imgs[2].height)),(255,255,255))
    resample=Image.BICUBIC

    if imgs[1].height == imgs[2].height:
        _im1 = imgs[1]
        _im2 = imgs[2]
        _im0 = imgs[0]
    else: # if (((imgs[1].height < imgs[2].height))) or also do the same for if imgs[2].height < imgs[1].height
        #_im2 = imgs[2].resize((int(imgs[2].width * imgs[1].height / imgs[2].height), int(imgs[1].height * imgs[2].width / imgs[1].width)), resample=resample)
        #_im2 = imgs[2].resize((imgs[2].width, int(imgs[2].height * imgs[1].width / imgs[2].width)), resample=resample)
        #_im2 = imgs[2].resize((imgs[2].width, ), resample=resample)
        
        # for very horizontal laytout
        #_im2 = ImageOps.contain(imgs[2],(imgs[1].width,imgs[1].height))

        # for square layout
        _im2 = ImageOps.contain(imgs[2],(imgs[1].width,imgs[1].height+imgs[0].height))
        _im1 = imgs[1]
        _im0 = imgs[0]

    # for having them side by side by side
    #dst = Image.new('RGB', (_im0.width+_im1.width + _im2.width, _im0.height+_im2.height),(255,255,255))

    # or for having them top & bottom then fitsmap on right side (overall square layout)
    dst = Image.new('RGB', (_im0.width+_im2.width,max(_im0.height+_im1.height,_im2.height)),(255,255,255))
    dst.paste(_im0, (0, 0))
    #dst.paste(_im1, (_im0.width, 0))
    dst.paste(_im1, (0, _im0.height))
    dst.paste(_im2, (_im0.width,0))
    filename = f'{kirk_names[c]}_{jades_names[c]}.jpg'
    dst.save(f'/home/robbler/research/Galaxies/form_images/{filename}')
    print(f'Saved file: {filename}')
    (i.close() for i in imgs)
    c+=1
    