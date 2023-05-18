# https://github.com/ysBach/SNU_AOpython/blob/main/chaps/01-fits-basic.ipynb
"""from IPython.core.interactiveshell import InteractiveShell
from IPython import get_ipython
#%config InlineBackend.figure_format = 'retina'
InteractiveShell.ast_node_interactivity = 'last_expr'
ipython = get_ipython()"""

from pathlib import Path
import numpy as np

from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.nddata import CCDData

from matplotlib import pyplot as plt
from matplotlib import rcParams

import ysfitsutilpy as yfu
import ysphotutilpy as ypu

plt.style.use('default')
rcParams.update({'font.size':12})

import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)

import _tool_visualization as vis

### Change here!!!
DATAPATH = Path('./Data/Cal_Data')
TMPDIR = Path('tmp')
TMPDIR.mkdir(exist_ok=True)

# In 3
### Change here !!!
allfits = list(DATAPATH.glob("M13*.fit")) 
allfits.sort()
print(f"Found {len(allfits)} fits files in {DATAPATH}.")
for _f in allfits:
    print(_f)

# In 4
for fpath in allfits:
    hdul = fits.open(fpath)
    hdul.info()

"""# In 5
print("Type of hdul:")
print(type(hdul), end='\n\n')

print("hdul.info():")
hdul.info()

# In 6
print("First 2x2 data of the 0-th element of hdul: hdul[0].data")
print(hdul[0].data[:2, :2], end='\n\n')
print("First 2x2 data of the 1-th element of hdul: hdul[1].data")
print(hdul[1].data[:2, :2], end='\n\n')

# In 7
print("First 10 header keywords for hdul[0].header")
print(list(hdul[0].header.keys())[:10], end='\n\n')

print("DATE-OBS (UT time of the **start** of the exposure) of hdul[0]: hdul[0].header['DATE-OBS']")
print(hdul[0].header['DATE-OBS'])
"""
# Close the HDU
hdul.close()

# In 8
# Only one file, wrap into loop
for fit in allfits:
    hdul = fits.open(allfits[0])

    # In 9
    obst = Time(hdul[0].header["DATE-OBS"])
    expt = hdul[0].header["EXPTIME"] * u.s
    print("The start of the observation time is  :", obst)
    print("The middle of the observation time is :", obst + expt/2)
    print("The end of the observation time is    :", obst + expt)



    fig, axs = plt.subplots(1, 1, figsize=(4, 3), sharex=False, sharey=False, gridspec_kw=None)
    vis.norm_imshow(axs, hdul[0].data, zscale=True)
    plt.tight_layout()

    hdul = fits.open(allfits[0])
    data = hdul[0].data
    newdata = np.sqrt(data[300:400, 650:750])
    print(newdata.shape)

    fig, axs = plt.subplots(1, 2, figsize=(5, 3), sharex=False, sharey=False, gridspec_kw=None)
    vis.norm_imshow(axs[0], newdata, zscale=True)
    vis.norm_imshow(axs[1], newdata, zscale=False, stretch='sqrt')

    plt.tight_layout()

    newhdu = fits.PrimaryHDU(data=newdata, header=hdul[0].header)
    newhdu.writeto(Path("tmp") / "test.fits", overwrite=True, output_verify='fix')

from astropy.nddata import Cutout2D
from astropy.wcs import WCS


def cut_ccd(ccd, position, size, mode="trim", fill_value=np.nan, warnings=True):
    """ Converts the Cutout2D object to proper CCDData."""
    cutout = Cutout2D(
        data=ccd.data, position=position, size=size,
        wcs=getattr(ccd, "wcs", WCS(ccd.header)),
        mode=mode, fill_value=fill_value, copy=True,
    )
    # Copy True just to avoid any contamination to the original ccd.

    nccd = CCDData(data=cutout.data, header=ccd.header.copy(),
                   wcs=cutout.wcs, unit=ccd.unit)
    ny, nx = nccd.data.shape
    nccd.header["NAXIS1"] = nx
    nccd.header["NAXIS2"] = ny
    nccd.header["LTV1"] = nccd.header.get("LTV1", 0) - cutout.origin_original[0]
    nccd.header["LTV2"] = nccd.header.get("LTV2", 0) - cutout.origin_original[1]

    return nccd

for fit in allfits:
    ccd = CCDData.read(fit, unit= "adu") # unit header needed for not example files
    cut = cut_ccd(ccd, position=(700, 350), size=(100, 100))
    print(cut.shape)
    cut.data = np.sqrt(cut.data)

    fig, axs = plt.subplots(1, 2, figsize=(5, 3), sharex=False, sharey=False, gridspec_kw=None)
    vis.norm_imshow(axs[0], cut.data, zscale=True)
    vis.norm_imshow(axs[1], cut.data, zscale=False, stretch='sqrt')
    plt.tight_layout()

    cut.write(Path("tmp") / "test.fits", overwrite=True, output_verify='fix')

    ccd = CCDData.read(allfits[0], unit= "adu")
    cut = yfu.imslice(ccd, trimsec="[650:749, 300:399]")
    print(cut.shape)
    cut.data = np.sqrt(cut.data)

    fig, axs = plt.subplots(1, 2, figsize=(5, 3), sharex=False, sharey=False, gridspec_kw=None)
    vis.norm_imshow(axs[0], cut.data, zscale=True)
    vis.norm_imshow(axs[1], cut.data, zscale=False, stretch='sqrt')
    plt.tight_layout()

print ("CCD over")

import astro_ndslice as nds
arr2d = np.arange(100).reshape(10,10)
a_slice = nds.slice_from_string('[1:3, 1:4]', fits_convention=True)
arr2d[a_slice]

a_slice = nds.slice_from_string('[1:10:2, 1:4]', fits_convention=True)
arr2d[a_slice]