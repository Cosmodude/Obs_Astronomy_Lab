from pathlib import Path
import astroalign as aa
from astropy.nddata import CCDData
from astropy.io import fits

DATAPATH = Path('./Data/Raw_Data')
TMPDIR = Path('tmp')
TMPDIR.mkdir(exist_ok=True)

### Astroalign 
### for every filter allign images one by one 
#### Throws errors : Big-endian buffer for ccd and reference stars less than min for data
def align_and_stack(array):
    print(len(array))
    stacked_image = array[0]
    for i in range(1, len(array)):
        print(i)
        target = array[0]#.newbyteorder()
        source = array[i]#.newbyteorder()
        # here comes the error, only with calibrated data
        aligned_source, footprint = aa.register( source, target , propagate_mask=True,detection_sigma=3) 
        # stack images
        stacked_image += aligned_source
    print(stacked_image)

    # Saving file
    DATAPATH = Path('./Data/Aligned_Raw/')
    fits.writeto(DATAPATH / (f"M13-{filter}" + '_aligned.fit'), stacked_image, header_obj[0], overwrite = True)
    return stacked_image

filter= "V"
allfits = list(DATAPATH.glob(f'M13*{filter}*.fit'))
allfits.sort()
data = []
ccd_array= []
header_obj = []
for idx, fit in enumerate(allfits):
    ccd_array.append(CCDData.read(fit,unit= "adu"))
    d, h = fits.getdata(fit, header = True)
    header_obj.append(h)
    data.append(d)
#print(ccd_array[0].wcs)

align_and_stack(ccd_array)


