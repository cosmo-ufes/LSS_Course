#%% IMPORT
import os
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import fits 
from astropy.table import Table
from astropy import units as u
from astropy.cosmology import Planck18
from astropy.coordinates import Distance as Distance
import pickle


#%% CATÁLOGO Y SELECCIÓN DE DATOS

os.chdir('/home/ftoscano/Doctorado/Quasars/')

with fits.open('DR16Q_v4.fits') as catalog:
    tabla_datos = Table(catalog[1].data)
    data = catalog[1].data

D_muestra = data[(data['Z']>=0.8) & (data['Z']<=2.2)]

#%% COORDENADAS


# Coordenadas con plt
plt.figure()
plt.subplot(111, projection = 'aitoff')
plt.grid(True)
plt.scatter((D_muestra['RA']-180.)*np.pi/180., (D_muestra['DEC'])*np.pi/180., s = 0.1,color='blue')
plt.title('QSOs sky position')
plt.show()

# Coordenadas con hp
nside = 2048
tam = hp.nside2npix(nside)
m = np.ones(tam)
hp.mollview(m, title='QSOs sky position',cbar=False, cmap='inferno')
hp.projscatter(D_muestra['RA'],D_muestra['DEC'],lonlat=True, s=0.1, color='orange',label='0, 360')
hp.graticule()
plt.legend()
plt.show()

hp.mollview(m, title='QSOs sky position',cbar=False, cmap='inferno')
hp.projscatter(D_muestra['RA']-180,D_muestra['DEC'],lonlat=True, s=0.1,color='green', label='-180,180')
hp.graticule()
plt.legend()
plt.show()

hp.mollview(m, title='QSOs sky position',cbar=False,cmap='inferno')
hp.projscatter(D_muestra['RA']-180,D_muestra['DEC'],lonlat=True, s=0.2,color='green', alpha=0.5,label='-180, 180')
hp.projscatter(D_muestra['RA'],D_muestra['DEC'],lonlat=True, s=0.2, color='orange', alpha = 0.5, label='0, 360')
hp.graticule()
plt.legend()
plt.show()

# Máscara de la Milky Way

#https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Lensing

os.chdir('/home/ftoscano/Doctorado/Mapas_CMB/')
maskalm_fname='rotated_alm._mask.dat'

with open(maskalm_fname, "rb" ) as f:  
   mask_alm=pickle.load( f)

m[mask_alm<0.75] = hp.pixelfunc.UNSEEN
hp.mollview(m, title='QSOs sky position with Milky Way mask',cbar=False, cmap='inferno')
hp.projscatter(D_muestra['RA']-180,D_muestra['DEC'],lonlat=True, s=0.2,color='green', alpha=0.5, label='-180, 180')
hp.projscatter(D_muestra['RA'],D_muestra['DEC'],lonlat=True, s=0.2, color='orange', alpha = 0.5, label='0, 360')
hp.graticule()
plt.legend()
plt.show()

hp.mollview(m, title='QSOs sky position with Milky Way mask',cbar=False, cmap='inferno')
hp.projscatter(D_muestra['RA'],D_muestra['DEC'],lonlat=True, s=0.2, color='orange', alpha = 0.5,label='0, 360')
hp.graticule()
plt.legend()
plt.show()

#%% HISTOGRAMAS

plt.hist(D_muestra['Z'], bins = 100, color = 'blue', alpha = 1., histtype='step',label = 'QSOs Redshift Distribution')
plt.legend()
plt.xlabel('z')
plt.ylabel('N(z)')
plt.show()

plt.hist(-D_muestra['M_I'], bins = 100, color = 'violet', alpha = 1., histtype='step',label = 'QSOs Absolute Magnitude Distribution')
plt.legend()
plt.xlabel('-M_I')
plt.ylabel('N(M_I)')
plt.show()

#Cortes en Magnitud

D_final = D_muestra[-D_muestra['M_I']>21]

plt.hist(-D_final['M_I'], bins = 100, color = 'violet', alpha = 1., histtype='step',label = 'QSOs Absolute Magnitude Distribution')
plt.legend()
plt.xlabel('-M_I')
plt.ylabel('N(M_I)')
plt.show()

#%% M_I vs z

plt.scatter(D_final['Z'],-D_final['M_I'], s=0.09, c='pink', marker="2",label='-M_I vs z')
plt.legend()
plt.xlabel('z')
plt.ylabel('-M_I')
plt.show()


#Subsamples

def f(a,x):
    return a+1.8*(x-0.8)

D_1 = D_final[-D_final['M_I']>f(25,D_final['Z'])]
#D_2 = D_final[-D_final['M_I']<f(22.3,D_final['Z'])]
D_2 = D_final[(-D_final['M_I']<f(22.3,D_final['Z'])) & ((-D_final['M_I']>f(21.3,D_final['Z'])))]
z_range = np.linspace(0.8,2.2,50)

plt.scatter(D_1['Z'],-D_1['M_I'], s=0.09, c='blue', marker="2",label='Greater -M_I vs z')
plt.scatter(D_2['Z'],-D_2['M_I'], s=0.09, c='magenta', marker="2",label='Lower -M_I vs z')
plt.plot(z_range,f(25,z_range), lw=1, c='black',linestyle='dashed',label='Greater fit')
plt.plot(z_range,f(22.3,z_range), lw=1, c='black',linestyle='dotted',label='Lower fit')
plt.legend()
plt.xlabel('z')
plt.ylabel('-M_I')
plt.show()




