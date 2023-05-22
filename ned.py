from astroquery.ned import Ned
import astropy.units as u
from astropy import coordinates

import sys

def band_filter(band, qualifier):
    if band.find('GALEX') != -1:
        if qualifier.find('diameter') != -1:
            return 1
        else:
            return -1
    if band.find('SDSS') != -1:
        if band.find('PSF') != -1:
            return 1
        else:
            return -1
    if band.find('2MASS') != -1:
        if qualifier.find('aperture') != -1:
            return 1
        else:
            return -1
    if band.find('WISE') != -1:
        if qualifier.find('Profile-fit') != -1:
            return 1
        else:
            return -1

tp = sys.argv[1]
with open("../sed/"+tp+"_oiii/"+tp+"_oiii_list") as f:
    lines = f.readlines()
    lines = [line.split() for line in lines]
    names = [line[0] for line in lines]
    ras = [float(line[1]) for line in lines]
    decs = [float(line[2]) for line in lines]

tt = open("../sed/"+tp+"_oiii/"+tp+"_oiii_list","w+")
try:
    for i in range(len(ras)):
        co = coordinates.SkyCoord(ra=ras[i],dec=decs[i],unit=(u.deg,u.deg))
        try:
            result_table = Ned.query_region(co, radius=0.0003*u.deg)
            result_table.sort('Photometry Points')
            pho = Ned.get_table(result_table['Object Name'][len(result_table)-1], table='photometry')
        except:
            print(names[i],ras[i],decs[i],sep=' ',file=tt)
            continue
        else:
            sed = pho[pho['NED Units']==b'Jy']['Observed Passband','Frequency',\
                'Flux Density','Upper limit of uncertainty','Lower limit of uncertainty',\
                'Qualifiers','Comments']
            out_f = open("../sed/"+tp+"_oiii/"+names[i]+".txt","w+")
            for j in range(len(sed)):
                if band_filter(str(sed[j]['Observed Passband']).replace(' ',''),\
                    str(sed[j]['Qualifiers']).replace(' ','')) == 1:
                    print(sed[j]['Frequency'],sed[j]['Flux Density'],\
                        sed[j]['Upper limit of uncertainty'],\
                            sed[j]['Lower limit of uncertainty'],file=out_f)
            out_f.close()
    tt.close()
except:
    for j in range(i,len(ras)):
        print(names[j],ras[j],decs[j],sep=' ',file=tt)
    tt.close()