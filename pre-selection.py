from astropy.io import fits
import pandas as pd
import numpy as np

cat = fits.open("../SDSS_dr14_properties.fits")
hdu1 = pd.DataFrame(cat[1].data).apply(
    pd.to_numeric, errors='ignore', downcast='float')

temp = hdu1[(hdu1['EWOIII5007'] / hdu1['e_EWOIII5007'] > 3)
            & (hdu1['EWHb-BR'] / hdu1['e_EWHb-BR'] > 3)
            & (hdu1['EWHb-NA'] / hdu1['e_EWHb-NA'] > 3)
            & (hdu1['EWHa-BR'] / hdu1['e_EWHa-BR'] > 3)
            & (hdu1['EWHa-NA'] / hdu1['e_EWHa-NA'] > 3)]

temp.to_csv(temp.csv)