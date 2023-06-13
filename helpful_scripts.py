from pathlib import Path
import astroalign as aa
from astropy.nddata import CCDData
from astropy.io import fits

DATAPATH = Path('./Data/Raw_Data')
TMPDIR = Path('tmp')
TMPDIR.mkdir(exist_ok=True)

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
#print(ccd_array)

### Astroalign 
### for every filter allign images one by one 
#### Throws errors : Big-endian buffer for ccd and reference stars less than min for data
def astrosync(array):
    print(len(array)-1)
    for i in range(0, len(ccd_array)-1):
        print(i)
        target = array[i]#.newbyteorder()
        source = array[i+1]#.newbyteorder()
        # here is the error only with calibrated data
        aligned_image, footprint = aa.register( source, target , propagate_mask=True,detection_sigma=3) 
        # set next target to alligned image
        array[i+1] = aligned_image
    print(aligned_image)

    # Saving file
    DATAPATH = Path('./Data/Aligned_Raw/')
    fits.writeto(DATAPATH / (f"M13-{filter}" + '_aligned.fit'), aligned_image, header_obj[0], overwrite = True)
    return aligned_image
astrosync(ccd_array)


