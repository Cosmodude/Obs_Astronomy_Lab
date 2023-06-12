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

from astropy import units as u
from astropy.nddata import CCDData, Cutout2D
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.stats import sigma_clipped_stats

from astroquery.mast import Catalogs
from astroquery.jplhorizons import Horizons

import ysfitsutilpy as yfu
import ysphotutilpy as ypu

from matplotlib import pyplot as plt
from matplotlib import rcParams
plt.style.use('default')
rcParams.update({'font.size':12})

from photutils.aperture import (CircularAperture, CircularAnnulus)
from photutils.detection import DAOStarFinder
from photutils.psf.groupstars import DAOGroup

from scipy.optimize import minimize

import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)

import _tool_visualization as vis

DATAPATH = Path('./Data/Cal_Data')
TMPDIR = Path('tmp')
TMPDIR.mkdir(exist_ok=True)

allfits = list(DATAPATH.glob("M13*V*.fit"))
allfits.sort()

# Only one file, wrap into loop

ccd = CCDData.read(allfits[0], unit= "adu")

# 1.1
for key in ["DATE-OBS", "EXPTIME", "FILTER"]:
    print(f"{key:10s} {str(ccd.header[key]):20s} {ccd.header.comments[key]}")

# 1.2
"""
objname = "4179"
observat = "B31"
t_obs = Time(ccd.header["DATE-OBS"]) + ccd.header["EXPTIME"] * u.s / 2
obj = Horizons(id=objname, location=observat, epochs=t_obs.jd)
q_obj = obj.ephemerides()


pos_sky = SkyCoord(q_obj["RA"][0], q_obj["DEC"][0], unit='deg')
pos_pix = pos_sky.to_pixel(wcs=ccd.wcs)
pos_pix = np.array([pos_pix[0], pos_pix[1]]).T

ap0 = CircularAperture(pos_pix, r=15)
an0 = CircularAnnulus(pos_pix, r_in=25, r_out=40)

fig, axs = plt.subplots(1, 1, figsize=(8, 5), sharex=False, sharey=False, gridspec_kw=None)

vis.norm_imshow(axs, ccd, zscale=True)
ap0.plot(color='r', lw=1, alpha=1, ax=axs)
an0.plot(color='w', lw=1, alpha=1, ax=axs)
plt.tight_layout()
plt.show();

phot_targ = ypu.apphot_annulus(ccd, ap0, an0)
m_targ = phot_targ['mag'][0]
dm_targ = phot_targ['merr'][0]
print(f"Instrumental magntiude = {m_targ:.3f} ± {dm_targ:.3f}")
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
