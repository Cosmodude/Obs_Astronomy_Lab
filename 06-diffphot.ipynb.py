"""
# %matplotlib notebook
from IPython.core.interactiveshell import InteractiveShell
from IPython import get_ipython
%config InlineBackend.figure_format = 'retina'
InteractiveShell.ast_node_interactivity = 'last_expr'
ipython = get_ipython()
"""

from pathlib import Path
import numpy as np

from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.nddata import CCDData, Cutout2D
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord

from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
from astroquery.jplhorizons import Horizons 
Gaia.ROW_LIMIT = -1
from matplotlib import pyplot as plt
from matplotlib import rcParams
plt.style.use('default')
rcParams.update({'font.size':12})

from photutils.centroids import centroid_com

import ysfitsutilpy as yfu
import ysphotutilpy as ypu

import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)

import _tool_visualization as vis
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

from astropy.table import Table
from photutils.psf.groupstars import DAOGroup

### Set yours here
DATAPATH = Path('./Data/Aligned_Cal')
TMPDIR = Path('tmp')
TMPDIR.mkdir(exist_ok=True)

### Change filters here manually
allfits = list(DATAPATH.glob("M13*V*.fit"))
allfits.sort()

###  Only one file, wrap into loop

ccd = CCDData.read(allfits[0], unit= "adu")

###  Cut outer regions of the image, keep only center
def cut_ccd(ccd, position, size, mode="trim", fill_value=np.nan, warnings=True):
    """ Converts the Cutout2D object to proper CCDData """
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

cutccd = cut_ccd(ccd,position=(2000,2000), size=(2000,2000))

### Query stars using DAO starfinder using circular aperture

from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture

avg, med, std = sigma_clipped_stats(cutccd.data)  # by default, 3-sigma 5-iteration.
finder = DAOStarFinder(threshold=4*std, fwhm=6.7, roundlo=-.5, roundhi=.5, exclude_border=True) 
sources = finder(cutccd.data - med) 

for col in sources.colnames:  
    sources[col].info.format = "%d" if col in ('id', 'npix') else '%.2f'
sources.pprint(max_width=76)  # astropy Table

fig, axs = plt.subplots(1, 1, figsize=(8, 5), sharex=False, sharey=False, gridspec_kw=None)

vis.norm_imshow(axs, cutccd.data, zscale=True)
pos = np.transpose((sources['xcentroid'], sources['ycentroid']))
aps = CircularAperture(pos, r=20.)
aps.plot(color='blue', lw=1, alpha=0.5)

plt.tight_layout()
plt.show();

"""
# 1.3 

pix_scale = 0.4*u.arcsec
center_xy = np.array(ccd.shape)/2
center_radec = ccd.wcs.wcs_pix2world(*center_xy, 0)
center_coo = SkyCoord(*center_radec, unit='deg')
fov_radius = np.sqrt((np.array(ccd.shape)**2).sum())/2 * pix_scale
q_ps = Catalogs.query_region(center_coo, radius=fov_radius, catalog="Panstarrs", 
                             data_release="dr2", table="mean")
# Change some column names for convenience.
q_ps["raMean"].name = "ra"
q_ps["decMean"].name = "dec"
q_ps["gMeanPSFMag"].name = "g"
q_ps["rMeanPSFMag"].name = "r"

# drop stars with unknown magnitudes
q_ps = q_ps.to_pandas().dropna(subset=["g", "r"])

# Calculate V and R, and their errors
q_ps["V"] = 0.006 + 0.474*q_ps["g"] + 0.526*q_ps["r"]
q_ps["R"] = -0.138 - 0.131*q_ps["g"] + 1.131*q_ps["r"]
q_ps["dV"] = np.sqrt(
    0.474**2*q_ps["gMeanPSFMagErr"]**2 
    + 0.526**2*q_ps["rMeanPSFMagErr"]**2 + 0.012**2
)
q_ps["dR"] = np.sqrt(
    0.131**2*q_ps["gMeanPSFMagErr"]**2
    + 1.131**2*q_ps["rMeanPSFMagErr"]**2 + 0.015**2
)
q_ps["dgr"] = np.sqrt(q_ps["gMeanPSFMagErr"]**2 + q_ps["rMeanPSFMagErr"]**2)

# Select only important columns
q2 = q_ps[["ra", "dec", "g", "r", "dgr", "V", "R", "dV", "dR"]].copy().reset_index(drop=True)

# Select only brighter than 22 mag
q2 = q2[(q2["V"] < 22) & (q2["R"] < 22)].copy()

# Calculate x, y position
coo = SkyCoord(q2["ra"], q2["dec"], unit='deg')
q2["x"], q2["y"] = ccd.wcs.wcs_world2pix(coo.ra, coo.dec, 0)

# Remove stars outside the image
q2 = q2[(q2["x"] > 20) & (q2["x"] < ccd.shape[1]-20) 
        & (q2["y"] > 20) & (q2["y"] < ccd.shape[0]-20)]

q2 = q2.reset_index(drop=True)

print(f"Total {len(q2)} stars from PS1 DR2")

avg, med, std = sigma_clipped_stats(ccd.data)  # by default, 3-sigma 5-iteration.
finder = DAOStarFinder(
    fwhm=4,  # In reality, FWHM must be measured a priori using, e.g., `ginga`
    threshold=4 * std, 
    sharplo=0.2, sharphi=1.0,   # default values 0.2 and 1.0
    roundlo=-1.0, roundhi=1.0,  # default values -1 and +1
    sigma_radius=1.5,           # default values 1.5
    ratio=1.0,                  # 1.0: circular gaussian
    exclude_border=True         # To exclude sources near edges
)

# The DAOStarFinder object ``finder`` gets at least one input: the image.
# Then it returns the astropy table which contains the aperture photometry results:
dao = finder(ccd.data - med)  
print(f"Total {len(dao)} stars from DAOStarFinder with threshold = {finder.threshold:.1f}")

# 1.3.1
"""