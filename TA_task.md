### https://github.com/Cosmodude/Obs_Astronomy_Lab/blob/main/TA_task.md
# Set 1

### 1 
https://docs.astropy.org/en/stable/io/fits/index.html
### 2
I did

### 3 
WCS - World Coordinate System;

I did

# Set 2 
Rules: 
    Never use for or while loop.

    You should not import any other packages.

    Answer to each problem must be a one-line of python code.

    For each problem, I gave hints. It is also homework for you to search for those on google.

```
# https://github.com/Cosmodude/Obs_Astronomy_Lab/blob/main/TA_task.py
# Used https://docs.astropy.org/en/stable/io/fits/index.html
# Solved this: 
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

```