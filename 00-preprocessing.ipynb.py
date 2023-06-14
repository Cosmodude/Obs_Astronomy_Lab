# In 2

import numpy as np

from astropy.io import fits                       # to handle FITS file
from astropy.visualization import ZScaleInterval  # to display data in z-scale

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from pathlib import Path                          # it is useful for directory-related tasks

# In 3
### Change HERE !!!
DATAPATH = Path('./Data')                      # path of data directory
RAWPATH  = DATAPATH/'Raw_Data'              # directory for raw data
CALPATH  = DATAPATH/'Aligned_Cal' 
ALIGNEDRAWPATH = DATAPATH/'Aligned_Raw'

# In 4 
# Change between: B,I,R,V manually in the following line's path and run script for each filter
list_obj = list(ALIGNEDRAWPATH.glob('M13-*V*.fit'))       # gathering raw fits files of targets
list_obj.sort()                                    # sorting
list_obj         # directory for saving the calibrated data

# In 5
hdul = fits.open(list_obj[0])                       # open fits file (HDUList)
hdul.info()

# In 6 
header      = hdul[0].header                        # get header from HDU
exptime_obj = header['EXPTIME']                     # exposure time (sec)
filter_obj  = header['FILTER']                      # filter
header

# In 7
data_obj = hdul[0].data                            # get 2-dimensional raw data from HDU            

# plot figure
fig = plt.figure(figsize = (8, 8)); fontsize = 15

ax  = fig.add_subplot()
vmin, vmax = ZScaleInterval().get_limits(data_obj) # to plot in z-scale
ax.imshow(data_obj, vmin = vmin, vmax = vmax)
print(list_obj[0])
ax.set_title(f'Raw Image ({str(list_obj[0]).split("/")[-1]})', fontsize = fontsize)  # There  was a typo : 1 not -1
plt.show()
print("First figure plotted!")

# In 8 
list_bias = list(RAWPATH.glob('calibration-*bias.fit'))
list_dark = list(RAWPATH.glob(f'calibration-*dk{int(exptime_obj)}.fit')) # we need the dark with same exposure time
list_flat = list(RAWPATH.glob(f'Flat-*{filter_obj}.fit'))                # we need the flat with same filter

# In 9 
list_bias.sort()
list_bias

# In 10 
list_dark.sort()
list_dark

# In 11 
list_flat.sort()
list_flat

# In 12 
data_bias = fits.getdata(list_bias[0])
data_dark = fits.getdata(list_dark[0])
data_flat = fits.getdata(list_flat[0])

dict_data = {'Bias'   : data_bias,
             'Dark'   : data_dark,
             'Flat'   : data_flat,
             'Object' : data_obj}

fig = plt.figure(figsize = (10, 10)); fontsize = 15
gridspec = GridSpec(nrows = 2, ncols = 2)

for i, keys in enumerate(dict_data.keys()):
    
    data = dict_data[keys]
    
    ax = fig.add_subplot(gridspec[i])
    vmin, vmax = ZScaleInterval().get_limits(data)
    ax.imshow(data, vmin = vmin, vmax = vmax)
    ax.axis('off')
    ax.set_title(keys, fontsize = fontsize)
    
plt.tight_layout()
plt.show()
print("Second figure plotted!")

# In 13 
array_bias = np.zeros((len(list_bias), np.shape(data_bias)[0], np.shape(data_bias)[1])) # Ndata x 4096 x 4096

for i, bias in enumerate(list_bias):
    data_bias = fits.getdata(bias)
    array_bias[i] = data_bias

data_mbias = np.median(array_bias, axis = 0) # master bias (median combining)

# In 14 
array_dark = np.zeros((len(list_dark), np.shape(data_dark)[0], np.shape(data_dark)[1]))

for i, dark in enumerate(list_dark):
    data_dark = fits.getdata(dark)
    array_dark[i] = data_dark

data_mdark = np.median(array_dark, axis = 0) - data_mbias # master dark
data_mdark[data_mdark<0] = 0.                             # Correct the negative values

# In 15 
array_flat = np.zeros((len(list_flat), np.shape(data_flat)[0], np.shape(data_flat)[1]))

# extract exposure time from flat images, and repeat 3.2.
exptime_flat = fits.getheader(list_flat[0])['EXPTIME']
list_fdark = list(RAWPATH.glob(f'calibration-*dk{int(exptime_flat)}.fit'))
list_fdark.sort()
data_fdark = fits.getdata(list_fdark[0])

array_fdark = np.zeros((len(list_fdark), np.shape(data_fdark)[0], np.shape(data_fdark)[1]))
for i, fdark in enumerate(list_fdark):
    data_fdark = fits.getdata(fdark)
    array_fdark[i] = data_fdark

data_mfdark = np.median(array_fdark, axis = 0) - data_mbias
data_mfdark[data_mfdark<0] = 0.                             

# master flat
array_flat = np.zeros((len(list_flat), np.shape(data_flat)[0], np.shape(data_flat)[1]))

for i, flat in enumerate(list_flat):
    data_flat = fits.getdata(flat)
    array_flat[i] = data_flat

data_mflat = np.median(array_flat, axis = 0) - data_mfdark - data_mbias # master flat
data_mflat /= np.max(data_mflat)                                        # normalization

# In 16
data_cal_obj = (data_obj - data_mbias - data_mdark) / data_mflat # calibrated data

# In 17
dict_data = {'Before calibration': data_obj,
             'After calibration' : data_cal_obj}

fig = plt.figure(figsize = (10, 10)); fontsize = 15
gridspec = GridSpec(nrows = 1, ncols = 2)

for i, keys in enumerate(dict_data.keys()):
    
    data = dict_data[keys]
    
    ax = fig.add_subplot(gridspec[i])
    vmin, vmax = ZScaleInterval().get_limits(data)
    ax.imshow(data, vmin = vmin, vmax = vmax)
    ax.axis('off')
    ax.set_title(keys, fontsize = fontsize)
    
plt.tight_layout()
plt.show()
print("Third figure plotted!")
# In 18 
for obj in list_obj:
    
    # image calibration for all images
    data_obj, header_obj = fits.getdata(obj, header = True)
    data_cal_obj = (data_obj - data_mbias - data_mdark) / data_mflat
    data_cal_obj -= data_cal_obj.min()
    print("File:", obj)
    print("Min data is", data_cal_obj.min())
    
    CALPATH.mkdir(exist_ok = True)
    fits.writeto(CALPATH / (str(obj).split('.')[0].split("\\")[-1] + '_cal.fit'), data_cal_obj, header_obj, overwrite = True)
    print("File saved!")