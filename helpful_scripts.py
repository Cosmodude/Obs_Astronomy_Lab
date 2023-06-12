from pathlib import Path
import astroalign as aa
from astropy.nddata import CCDData
from astropy.io import fits

DATAPATH = Path('./Data/Cal_Data')
TMPDIR = Path('tmp')
TMPDIR.mkdir(exist_ok=True)

allfits = list(DATAPATH.glob("M13*R*.fit"))
allfits.sort()
data = []
for fit in allfits:
    data.append(fits.getdata(fit))
ccd_array= []
for fit in allfits:
    ccd_array.append(CCDData.read(fit,unit= "adu"))
#print(ccd_array)

### Astroalign 
### for every filter allign images one by one 
#### Throws errors : Big-endian buffer for ccd and reference stars less than min for data
def astrosync():
    for i in range(0, len(ccd_array)-1):
        print(i)
        target = ccd_array[0]#.newbyteorder()
        source = ccd_array[1]#.newbyteorder()
        # hear is the error
        aligned_images_arr, footprint = aa.register( source, target , propagate_mask=True,detection_sigma=2) 
    print(aligned_images_arr)



