# https://github.com/ysBach/SNU_AOpython/blob/main/chaps/02-center.ipynb
"""from IPython.core.interactiveshell import InteractiveShell
from IPython import get_ipython
%config InlineBackend.figure_format = 'retina'
InteractiveShell.ast_node_interactivity = 'last_expr'
ipython = get_ipython()"""

# In 30
from pathlib import Path
import numpy as np

from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.nddata import CCDData, Cutout2D
from astropy.stats import sigma_clipped_stats

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

DATAPATH = Path('./Data/Cal_Data')
TMPDIR = Path('tmp')
TMPDIR.mkdir(exist_ok=True)

# In 38
allfits = list(DATAPATH.glob("M13*.fit"))
allfits.sort()

# Only for one file, wrap into loop
for fit in allfits:
    ccd = CCDData.read(fit, unit= "adu")
    cut = Cutout2D(ccd, position=(270,320), size=(95,100))

    fig, axs = plt.subplots(1, 1, figsize=(3, 4), sharex=False, sharey=False, gridspec_kw=None)
    vis.norm_imshow(axs, cut.data, zscale=True)
    plt.tight_layout()

    # In 40
    avg, med, std = sigma_clipped_stats(cut.data)  # by default, 3-sigma 5-iteration.
    thresh_3sig = med + 3 * std
    mask_3sig = (cut.data < thresh_3sig)
    center = centroid_com(data=cut.data, mask=mask_3sig)
    print(center)
    print(cut.to_original_position(center))

    fig, axs = plt.subplots(1, 1, figsize=(3, 4), sharex=False, sharey=False, gridspec_kw=None)
    vis.norm_imshow(axs, mask_3sig.astype(int))
    vis.norm_imshow(axs, cut.data, alpha=0.4, zscale=True)
    axs.plot(*center, 'rx')
    plt.tight_layout()

    # In 41 
    from photutils.detection import DAOStarFinder
    from photutils.aperture import CircularAperture

    avg, med, std = sigma_clipped_stats(ccd.data)  # by default, 3-sigma 5-iteration.
    finder = DAOStarFinder(threshold=5.*std, fwhm=4, exclude_border=True) 
    sources = finder(ccd.data - med)  
    for col in sources.colnames:  
        sources[col].info.format = "%d" if col in ('id', 'npix') else '%.2f'
    sources.pprint(max_width=76)  # astropy Table

    fig, axs = plt.subplots(1, 1, figsize=(8, 5), sharex=False, sharey=False, gridspec_kw=None)

    vis.norm_imshow(axs, ccd.data, zscale=True)
    pos = np.transpose((sources['xcentroid'], sources['ycentroid']))
    aps = CircularAperture(pos, r=10.)
    aps.plot(color='blue', lw=1.5, alpha=0.5)

    plt.tight_layout()
    plt.show();

    # In 42
    sources["dist"] = np.sqrt(np.sum((pos - np.array([270, 320]))**2, axis=1))
    sources.sort("dist")
    sources[0]