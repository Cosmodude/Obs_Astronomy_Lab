from pathlib import Path
import astroalign as aa
from astropy.nddata import CCDData

DATAPATH = Path('./Data/Cal_Data')
TMPDIR = Path('tmp')
TMPDIR.mkdir(exist_ok=True)

allfits = list(DATAPATH.glob("M13*V*.fit"))
allfits.sort()
ccd= []
for fit in allfits:
    ccd.append(CCDData.read(fit))
print(ccd)

### Astroalign 
### for every filter allign them one by one 
def astrosync(images):
    aligned_image, footprint = aa.register(source_image, target_image)
    pass;
        