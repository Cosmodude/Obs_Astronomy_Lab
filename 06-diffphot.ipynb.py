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
#DATAPATH = Path('./Data/Cal_Data')
TMPDIR = Path('tmp')
TMPDIR.mkdir(exist_ok=True)

### Change filters here manually
allfits = list(DATAPATH.glob("M13*V*.fits"))
allfits.sort()

###  Only one file, wrap into loop

ccd = CCDData.read(allfits[0], unit = 'u')

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

### Find stars with DAOStarfinder using circular aperture

from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture

avg, med, std = sigma_clipped_stats(cutccd.data)  # by default, 3-sigma 5-iteration.
finder = DAOStarFinder(threshold=4*std, fwhm=6.7, roundlo=-.5, roundhi=.5, exclude_border=True) 
sources = finder(cutccd.data - med) 

for col in sources.colnames:  
    sources[col].info.format = "%d" if col in ('id', 'npix') else '%.2f'
#sources.pprint(max_width=76)  # astropy Table

fig, axs = plt.subplots(1, 1, figsize=(8, 5), sharex=False, sharey=False, gridspec_kw=None)

vis.norm_imshow(axs, cutccd.data, zscale=True)
pos = np.transpose((sources['xcentroid'], sources['ycentroid']))
aps = CircularAperture(pos, r=20.)
aps.plot(color='blue', lw=1, alpha=0.5)

plt.tight_layout()
plt.show();

# Find center of the image 

print(cutccd.wcs)
pix_scale = 0.4*u.arcsec

center_xy = np.array(cutccd.shape)/2
center_radec = cutccd.wcs.wcs_pix2world(*center_xy, 0)
center_coo = SkyCoord(*center_radec, unit='deg')

width, height = np.array(cutccd.shape)*pix_scale

print("\nCoordinate of the center of the image:\n", center_coo)

fov_radius = np.sqrt((np.array(cutccd.shape)**2).sum())/2 * pix_scale


# Query stars from the PanStarrs catalog

q_ps = Catalogs.query_region(center_coo, radius=fov_radius, catalog="Panstarrs", 
                             data_release="dr2", table="mean")
# Change some column names for convenience.
q_ps["raMean"].name = "ra"
q_ps["decMean"].name = "dec"
q_ps["gMeanPSFMag"].name = "g"
q_ps["rMeanPSFMag"].name = "r"
print("Length of queried catalog", len(q_ps))


# Dropping stars!!!
q_ps = q_ps.to_pandas().dropna(subset=["g", "r"])

# Calculate V and R magnitudes and B-V
q_ps["V"] = 0.006 + 0.474*q_ps["g"] + 0.526*q_ps["r"]
q_ps["R"] = -0.138 - 0.131*q_ps["g"] + 1.131*q_ps["r"]
q_ps["B-V"] = 0.207 + 1.113*(q_ps["g"]-q_ps["r"])
q_ps["dV"] = np.sqrt(
    0.474**2*q_ps["gMeanPSFMagErr"]**2 
    + 0.526**2*q_ps["rMeanPSFMagErr"]**2 + 0.012**2
)
q_ps["dR"] = np.sqrt(
    0.131**2*q_ps["gMeanPSFMagErr"]**2
    + 1.131**2*q_ps["rMeanPSFMagErr"]**2 + 0.015**2
)

q_ps["dgr"] = np.sqrt(q_ps["gMeanPSFMagErr"]**2 + q_ps["rMeanPSFMagErr"]**2)

# Select only necessary columns
q2 = q_ps[["ra", "dec", "g", "r", "dgr", "V", "R", "B-V", "dV", "dR"]].copy().reset_index(drop=True)

# Select only brighter than 22 mag
q2 = q2[(q2["V"] < 22) & (q2["R"] < 22)].copy().reset_index(drop=True)

# Calculate x, y position
coo = SkyCoord(q2["ra"], q2["dec"], unit='deg')
q2["x"], q2["y"] = cutccd.wcs.wcs_world2pix(coo.ra, coo.dec, 0)

# Remove stars outside the image
q2 = q2[(q2["x"] > 20) & (q2["x"] < cutccd.shape[1]-20) 
        & (q2["y"] > 20) & (q2["y"] < cutccd.shape[0]-20)]

q2 = q2.reset_index(drop = True)

# Grouping stars function 
def group_stars(table, crit_separation, xcol='x', ycol='y', index_only=True):
    if not isinstance(table, Table):
        table = Table.from_pandas(table)
    tab = table.copy()
    
    tab[xcol].name = "x_0"
    tab[ycol].name = "y_0"
    try:
        gtab = DAOGroup(crit_separation=crit_separation)(tab)
    except IndexError:
        gtab = tab
        gtab["group_id"] = []
        gtab["id"] = []
    if not index_only:
        gtab["x_0"].name = xcol
        gtab["y_0"].name = ycol
        return gtab
    else:
        gid, gnum = np.unique(gtab["group_id"], return_counts = True)
        gmask = gid[gnum != 1]
        grouped_rows = []
        for i, gid in enumerate(gtab["group_id"]):
            if gid in gmask:
                grouped_rows.append(i)
        return grouped_rows

# Dropping close stars 
rows2rm = group_stars(q2, crit_separation=10) 
q2_close = q2.drop(rows2rm, axis=0).reset_index(drop=True)

fig, axs = plt.subplots(1, 1, figsize=(8, 5), sharex=False, sharey=False, gridspec_kw=None)

aps_all = CircularAperture(np.array([q2["x"], q2["y"]]).T, r=15)
aps_close = CircularAperture(np.array([q2_close["x"], q2_close["y"]]).T, r=15)

vis.norm_imshow(axs, cutccd.data, zscale=True)

#Plotting all + close stars
aps_all.plot(color="r", lw=1, alpha=0.5, ax=axs)
aps_close.plot(color="blue", lw=1, alpha=0.5, ax=axs)

plt.title("Close stars")
print("Plotting 2nd plot")
plt.tight_layout()
plt.show();


rows2rm = group_stars(sources, crit_separation=10, xcol="xcentroid", ycol="ycentroid")
sources = sources.to_pandas().drop(rows2rm, axis=0).reset_index(drop=True)

# Plotting DAO + close stars
aps_dao = CircularAperture(np.array((sources['xcentroid'], sources['ycentroid'])).T, r=12)

fig, axs = plt.subplots(1, 1, figsize=(8, 5), sharex=False, sharey=False, gridspec_kw=None)

vis.norm_imshow(axs, cutccd.data, zscale=True)

aps_dao.plot(color="r", lw=1, alpha=0.5, ax=axs)
aps_close.plot(color="b", lw=1, alpha=0.5, ax=axs)
for c, l in zip("rb", ["DAO", "PS1"]):
    axs.plot([], [], c=c, label=l)
axs.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)

plt.tight_layout()
plt.title("DAO + PS1")
print("Plotting 3rd plot")
plt.show();

q2_close["xcentroid"], q2_close["ycentroid"] = None, None

for i, row in q2_close.iterrows():
    x, y = row["x"], row["y"]
    # Find closest DAOStarFinder source with distance no larger than 5 pix
    dx = sources["xcentroid"] - x
    dy = sources["ycentroid"] - y
    candidate = sources[(abs(dx) < 5) & (abs(dy) < 5)]
    if len(candidate) == 0:
        continue

    # Find the closest one
    if len(candidate) > 1:
        candidate = candidate[np.argmin(np.hypot(dx, dy))]
    
    q2_close.loc[i, "xcentroid"] = candidate["xcentroid"].iloc[0]
    q2_close.loc[i, "ycentroid"] = candidate["ycentroid"].iloc[0]
q2_close = q2_close.dropna(axis=0, subset="xcentroid").reset_index(drop=True)


### Working with aperture
from photutils.aperture import (CircularAperture, CircularAnnulus)
q2_xy = np.array([q2_close["xcentroid"], q2_close["ycentroid"]]).T
ap_stars = CircularAperture(q2_xy, r=10)
an_stars = CircularAnnulus(q2_xy, r_in=20, r_out=30)
phot = ypu.apphot_annulus(cutccd, ap_stars, an_stars).drop(["xcenter", "ycenter"], axis=1)
phot = q2_close.join([phot])

print(phot)

# Instrumental magnitude + B-V diagram
plt.plot(phot["B-V"], phot["mag"],"o", markersize=2)
plt.gca().invert_yaxis()
plt.xlabel("B-V")
plt.ylabel("instrumental V-mag")
plt.show

# Real magnitude +  B-V diagram
mcat, dmcat = phot["V"], phot["dV"]
mobs, dmobs = phot["mag"], phot["merr"]
color = phot["g"] - phot["r"]
dmtot = np.sqrt(dmcat**2 + dmobs**2)
mcat.fillna(0, inplace=True)
dmcat.fillna(0, inplace=True)
mobs.fillna(0, inplace=True)
dmobs.fillna(0, inplace=True)

# Standartization
zeropt = np.average(mobs - mcat, weights=1/dmtot**2)
mag = phot["mag"]-zeropt
color = phot["B-V"]

#print(mag)

plt.plot(color, mag, "o", markersize=2)
plt.gca().invert_yaxis()
plt.xlabel("B-V")
plt.ylabel("V-magnitude")
plt.show


### Fitting the curve
from scipy.optimize import curve_fit

def linf(x, a, b):
    return a + b*x

# === Calculate zero point and errors
zeropt = np.average(mobs - mcat, weights=1/dmtot**2)
dzeropt = np.max([1/np.sqrt(np.sum(1/dmtot**2)), np.std(mobs - mcat, ddof=1)/np.sqrt(len(mcat))])
dmtot2 = np.sqrt(dmtot**2 + dzeropt**2)

# === Find fitting lines
# Search for the usage of scipy.optimize.curve_fit.
poptm, _ = curve_fit(linf, mcat, mobs, sigma=dmobs, absolute_sigma=True)
poptc, _ = curve_fit(linf, color, mobs-mcat, sigma=dmtot2, absolute_sigma=True)

# === Plot
# --- Set some useful things
errkw = dict(marker="", ls="", ecolor="gray", elinewidth=0.5)

# --- Main plot with error bars and fitting lines
fig, axs = plt.subplot_mosaic("mc\nmc\nrc", figsize=(10, 5))
# m = magnitudes, c=colors, r=residuals
axs["m"].plot(mcat, mobs, "k.", ms=5)
axs["m"].errorbar(mcat, mobs, xerr=dmcat, yerr=dmobs, **errkw)
axs["r"].plot(mcat, mobs - mcat - zeropt, "k.", ms=5)
axs["r"].errorbar(mcat, mobs - mcat - zeropt, xerr=dmcat, yerr=dmtot2, **errkw)

mm = np.array(axs["m"].get_xlim())

# Fitted lines
axs["m"].plot(mm, zeropt + mm, "r-", lw=1, label=f"Z = {zeropt:+.3f} (fix slope)")
axs["m"].plot(mm, linf(mm, *poptm), "b:", lw=1, label=f"y = {poptm[1]:.3f}x {poptm[0]:+.3f}")

# --- Some codes to make the plot prettier
#axs["m"].axhline(m_targ, color="k", lw=1, label=f"Target v = {m_targ:.3f} ± {dm_targ:.3f}")
#axs["m"].hlines([m_targ+dm_targ, m_targ-dm_targ], *axs["m"].get_xlim(), color="k", lw=1, ls=":")

for i, row in phot.iterrows():
    axs["m"].text(row["V"], row["mag"], i, fontsize=8)
    axs["r"].text(row["V"], row["mag"] - row["V"] - zeropt, i, fontsize=8)
    

axs["m"].set(xlim=mm, ylabel="v (Instrumental V-mag)")
axs["r"].set(xlim=mm, ylim=np.array([-1, 1])*np.max(np.abs(axs["r"].get_ylim())), 
             ylabel="v - V_PS1 - Z", xlabel="PS1 V-mag")

axs["r"].axhline(0, color="k", lw=1)
axs["r"].hlines([dzeropt, -dzeropt], *mm, color="k", lw=1, ls=":" )

axs["m"].legend(fontsize=10)

plt.tight_layout()
plt.show()