# https://github.com/Cosmodude/Obs_Astronomy_Lab/blob/main/TA_task.py
# Used https://docs.astropy.org/en/stable/io/fits/index.html

### Task 1

import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.nddata import CCDData

np.random.seed(123)  # legacy function
SAVEPATH = Path("./Problem_products")  # <-- You may tune here for your computer
data = np.random.rand(100, 100) * 100 # Creates 100x100 array of random numbers
#print(data)

### Prob 1
hdu = fits.PrimaryHDU(data=data)
### Prob 2
hdu.data = hdu.data.astype(np.float32)
### Prob 3
hdu.writeto(SAVEPATH / "test.fits", overwrite=True)
### Prob 4
hdul = fits.open(SAVEPATH / "test.fits")
#print(hdul[0].data)

### Prob 5
print("Prob 5:")
hdul.info()  # prints result itself
### Prob 6 
print()
print("Prob 6:")
print(np.testing.assert_allclose(hdul[0].data,data)) # rises error if not equal, returns None
### Prob 7 
hdr_hdu= hdul[0].header
"""
print()
print("Prob 7:")
print(hdr_hdu)  # prints badly in terminal
"""
### Prob 8
ccd = CCDData.read(SAVEPATH / "test.fits",unit="parsec")
### Prob 9
np.testing.assert_allclose(ccd.data, data) # rises error if not equal, returns None
### Prob 10 
hdr_ccd = ccd.header
"""
print()
print("Prob 10:")
print(hdr_ccd) # prints badly in terminal
"""
### Prob 11
print()
print("Prob 11:")
print(hdr_ccd)
print(hdr_hdu)
# Prints same structure




### Task 2 Software & Tools

print()
print("Task 2")
import numpy
import scipy
import astropy
import pandas
import ccdproc
import photutils 
import specutils 

import astroscrappy
import matplotlib

print(ccdproc.__version__)
print(matplotlib.__version__)
print(numpy.__version__)
print(scipy.__version__)
print(astropy.__version__)
print(pandas.__version__)
print(photutils.__version__)
print(specutils.__version__)
print(astroscrappy.__version__)