from astropy.io import fits
import pandas as pd
import numpy as np

temp = pd.read_csv('temp.csv')

#fetching parameters 
name, ra, dec, redshift, plate, fiber, mjd, \
logl_oiii, err_logl_oiii, logl_hb_broad, \
err_logl_hb_broad, logl_hb_narrow, \
err_logl_hb_narrow, logl_ha_broad, \
err_logl_ha_broad, \
logl_ha_narrow, err_logl_ha_narrow, \
logl_nii, err_logl_nii, \
logl_mgii_broad, err_logl_mgii_broad, \
logl_mgii_narrow, err_logl_mgii_narrow, \
ew_oiii, err_ew_oiii, \
ew_hb_broad, err_ew_hb_broad, \
ew_mgii_broad, err_ew_mgii_broad, \
ew_mgii_narrow, err_ew_mgii_narrow, \
fwhm_hb_broad, err_fwhm_hb_broad, \
fwhm_mgii_broad, err_fwhm_mgii_broad, \
log_bh, err_log_bh, edd, log_5100, err_log_5100 \
= temp['SDSS'], temp['RAJ2000'], temp['DEJ2000'], \
temp['z'], temp['Plate'], temp['Fiber'], temp['MJD']
temp['logLOIII5007'], temp['e_logLOIII5007'], \
temp['logLHb-BR'], temp['e_logLHb-BR'], temp['logLHb-NA'], \
temp['e_logLHb-NA'], temp['logLHa-BR'], temp['e_logLHa-BR'], \
temp['logLHa-NA'], temp['e_logLHa-NA'], temp['logLNII6585'], \
temp['e_logLNII6585'], ['logLMgII-BR'], temp['e_logLMgII-BR'], \
temp['logLMgII-NA'], temp['e_logLMgII-NA'], temp['EWOIII5007'], \
temp['e_EWOIII5007'], temp['EWHb-BR'], temp['e_EWHb-BR'], temp['EWMgII-BR'], \
temp['e_EWMgII-BR'], temp['EWMgII-NA'], temp['e_EWMgII-NA'], temp['FWHMHb-BR'], \
temp['e_FWHMHb-BR'], temp['FWHMMgII-BR'], temp['e_FWHMMgII-BR'], temp['logMBH'], \
temp['e_logMBH'], temp['LogREdd'], temp['logL5100'], temp['e_logL5100']

log_bol, err_log_bol = [], []
for index, row in temp.iterrows():
    log_bol.append(row['logLbol'])
    if row['z'] < 0.8:
        err_log_bol.append(row['e_logL5100'])
    if 0.8 < row['z'] < 1.9:
        err_log_bol.append(row['e_logL3000'])
    if 1.9 < row['z']:
        err_log_bol.append(row['e_logL1350'])

# # oiii/hb_broad
# ratio = [oiii - hb for oiii,
#          hb in zip(logl_oiii, logl_hb_broad)]  
# ratio_error = [np.sqrt(oiii**2+hb**2) for oiii,
#          hb in zip(err_logl_oiii, err_logl_hb_broad)]  # error of ratio

# # hb_narrow/hb_broad
# ratio_hb = [n-b for n, b in zip(logl_hb_narrow, logl_hb_broad)]
# ratio_nii = [nii-ha for nii,
#              ha in zip(logl_nii, logl_ha_narrow)]  # nii/ha_narrow
# ratio_oiii = [oiii-hb for oiii,
#               hb in zip(logl_oiii, logl_hb_narrow)]  # oiii/hb_narrow
# ratio_balmer = [ha-hb for ha, hb in zip(logl_ha_broad, logl_hb_broad)]

# # the lowest 5% and the highest 5% of oiii/hb_broad
# # l, h are the two critical values
# l, h = np.percentile(ratio, 5, method='lower'), np.percentile(ratio, 95, method='higher')
# m = np.mean(ratio)

# # balmer decrement: select sources with balmer decrement in a narrow range
# mean_bmr = np.mean(ratio_balmer)
# a, b = mean_bmr-0.2, mean_bmr+0.2 # a, b are the two critical values
# c = np.mean(ratio_balmer)

# target_file = open("./low_oiii_list", "w+")
# for idx, r in enumerate(ratio):
#     if r <= l and ratio_balmer[idx] <= b and ratio_balmer[idx] >= a:
#         print(name[idx], ra[idx], dec[idx], plate[idx], fiber[idx],
#               mjd[idx], redshift[idx], r, sep=' ', file=target_file)
# target_file.close()

# target_file = open("./high_oiii_list", "w+")
# for idx, r in enumerate(ratio):
#     if r >= h and ratio_balmer[idx] <= b and ratio_balmer[idx] >= a:
#         print(name[idx], ra[idx], dec[idx], plate[idx], fiber[idx],
#               mjd[idx], redshift[idx], r, sep=' ', file=target_file)
# target_file.close()

# target_file = open("./norm_oiii_list", "w+")
# for idx, r in enumerate(ratio):
#     if r > l and r < h and ratio_balmer[idx] <= b and ratio_balmer[idx] >= a:
#         print(name[idx], ra[idx], dec[idx], plate[idx], fiber[idx],
#               mjd[idx], redshift[idx], r, sep=' ', file=target_file)
# target_file.close()

# target_file = open("./low_edd_list", "w+")
# for idx, r in enumerate(ratio):
#     if r >= h and ratio_balmer[idx] <= b and ratio_balmer[idx] >= a and edd[idx] <= -2:
#         print(name[idx], ra[idx], dec[idx], plate[idx], fiber[idx],
#               mjd[idx], redshift[idx], r, sep=' ', file=target_file)
# target_file.close()

# target_file = open("./low_bpt_list", "w+")
# for idx, r in enumerate(ratio):
#     if r <= l and ratio_balmer[idx] <= b and ratio_balmer[idx] >= a:
#         if logl_nii[idx] != 0 and (logl_oiii[idx]-logl_hb_narrow[idx]) \
#                 < (0.61 / ((logl_nii[idx]-logl_ha_narrow[idx]) - 0.05) + 1.3):
#             print(name[idx], ra[idx], dec[idx], plate[idx], fiber[idx],
#                   mjd[idx], redshift[idx], r, sep=' ', file=target_file)
# target_file.close()

# with open("./low_oiii_list") as f:
#     lines = f.readlines()
#     lines = [line.split() for line in lines]
#     low_oiii_cross_id = [line[3]+'-'+line[4]+'-'+line[5] for line in lines]

# with open("./high_oiii_list") as f:
#     lines = f.readlines()
#     lines = [line.split() for line in lines]
#     high_oiii_cross_id = [line[3]+'-'+line[4]+'-'+line[5] for line in lines]

# with open("./norm_oiii_list") as f:
#     lines = f.readlines()
#     lines = [line.split() for line in lines]
#     norm_oiii_cross_id = [line[3]+'-'+line[4]+'-'+line[5] for line in lines]

# with open("./low_edd_list") as f:
#     lines = f.readlines()
#     lines = [line.split() for line in lines]
#     low_edd_cross_id = [line[3]+'-'+line[4]+'-'+line[5] for line in lines]

# with open("./low_bpt_list") as f:
#     lines = f.readlines()
#     lines = [line.split() for line in lines]
#     low_bpt_cross_id = [line[3]+'-'+line[4]+'-'+line[5] for line in lines]

# index = [p+'-'+f+'-'+m for p,f,m in zip(plate,fiber,mjd)]
# data = {
#     'NAME':name,
#     'RA':ra,
#     'DEC':dec,
#     'REDSHIFT':redshift,

#     'PLATE':plate,
#     'FIBER':fiber,
#     'MJD':mjd,

#     'RATIO':ratio,
#     'RATIO_ERROR':ratio error,

#     'Lums_OIII':logl_oiii,
#     'Lums_HB_BROAD':logl_hb_broad,
#     'Lums_HB_NARROW':logl_hb_narrow,
#     'Lums_MGII_BROAD':logl_mgii_broad,

#     'EW_OIII':ew_oiii,
#     'EW_OIII_ERR':err_ew_oiii,
#     'EW_HB_BROAD':ew_hb_broad,
#     'EW_HB_BROAD_ERR':err_ew_hb_broad,
#     'EW_MGII_BROAD':ew_mgii_broad,
#     'EW_MGII_BROAD_ERR':err_ew_mgii_broad,

#     'EW_MGII_NARROW':ew_mgii_narrow,
#     'EW_MGII_NARROW_ERR':err_ew_mgii_narrow,

#     'FWHM_HB_BROAD':fwhm_hb_broad,
#     'FWHM_HB_BROAD_ERR': err_fwhm_hb_broad,
#     'FWHM_MGII_BROAD': err_fwhm_mgii_broad,
#     'FWHM_MGII_BROAD_ERR': err_fwhm_mgii_broad,

#     'MBH':log_bh,
#     'MBH_ERR':err_log_bh,

#     'Lbol':log_bol,
#     'Lbol_ERR':err_log_bol,

#     'Lums_5100':log_5100,
#     'Lums_5100_ERR':err_log_5100,

#     'Edd':edd,

#     'Lums_NII':logl_nii,
#     'Lums_HA_BROAD':logl_ha_broad,
#     'Lums_HA_NARROW':logl_ha_narrow,

#     'RATIO_HB':ratio_hb,
#     'RATIO_NII':ratio_nii,
#     'RATIO_OIII':ratio_oiii,
#     'RATIO_BALMER':ratio_balmer
# }
# catalog = pd.DataFrame(data, index)

# low_oiii_cat = catalog.loc[low_oiii_cross_id]
# high_oiii_cat = catalog.loc[high_oiii_cross_id]
# norm_oiii_cat = catalog.loc[norm_oiii_cross_id]
# low_edd_cat = catalog.loc[low_edd_cross_id]
# low_bpt_cat = catalog.loc[low_bpt_cross_id]

# print(len(high_oiii_cat),len(low_oiii_cat),len(norm_oiii_cat),sep=', ')
# print(len(low_edd_cat),len(low_bpt_cat),sep=', ')
