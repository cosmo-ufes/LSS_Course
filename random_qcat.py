# IMPORT MODULES/LIBRARIES __________________________________________________________________

import numpy as np # work with arrays, linear algebra, fourier transform, matrices, etc
import matplotlib.pyplot as plt # graphic facilities
from astropy.io import fits # read fits files
import healpy as hp
import pandas as pd




# HERE WE WRITE A SUBSET OF THE QSO CATALOGUE TO USE LESS RAM 
qso_data='data/DR16Q_v4.fits'

zmin=0.8
zmax=2.2

fileout='data/subsets/qso_subset_08_22.dat'

def qso_sub(qso_data,zmin,zmax,fileout):

    with fits.open(qso_data) as catalog:
        data0 = catalog[1].data

    data=data0[(data0['Z']>zmin)&(data0['Z']<zmax)] #---> z slice

    datasave=np.vstack([data['RA'],data['DEC'],data['Z']])                                                                                                              
    np.savetxt(fileout,datasave.T)

    print('qso number',len(data['RA']))

#qso_sub(qso_data,zmin,zmax,fileout) #call only once per sample

#___________________________________________________________________

# READ SUBSET OF QSO CATALOGUE______________________________________
colnames=['ra','dec','z']
D = pd.read_csv(fileout,delimiter=' ', names=colnames)

nside=32 #select healpy map resolution (pixel size)

background_healpix_map=np.zeros(hp.nside2npix(nside)) #create empty map to plot coordinates
hp.mollview(background_healpix_map,title='Angular distribution of QSO sample',cbar=False) #plot background map
hp.graticule()
hp.projscatter(D['ra'],D['dec'],lonlat=True,s=0.05,color='red') #plot RA DEC of QSO catalogue
plt.show()

plt.hist(D['z'],density=True,bins=25,histtype='step',label='QSO z distribution')
plt.legend() 
plt.xlabel('z')
plt.show()
# CREATE QSO ANGULAR MASK___________________________________________________

def make_ang_mask(nside,Dra,Ddec):

    mask=np.zeros(hp.nside2npix(nside))

    veclist=hp.ang2vec(Dra,Ddec,lonlat=True)

    pix=hp.vec2pix(nside,veclist[:,0],veclist[:,1],veclist[:,2])

    mask[pix]=1

    return(mask)

qso_ang_mask=make_ang_mask(nside,D['ra'],D['dec'])

hp.mollview(qso_ang_mask,title='Angular QSO mask')
plt.show()

hp.mollview(qso_ang_mask,title='Angular QSO mask and QSO positions')
hp.graticule()
hp.projscatter(D['ra'],D['dec'],lonlat=True,s=0.05,color='red') #plot RA DEC of QSO catalogue
plt.show()


#---------------------------------------------------------
# HOMOGENEUS RANDOM POSITIONS ON THE MASK AREA________________________________________________

#Nran=len(D['z'])
Nran=10000

RAran = np.zeros(Nran)
DECran = np.zeros(Nran)

nr=0

while(nr<Nran):

    RArand = np.random.uniform(0,360.,1)

    cosdec = np.random.uniform(-1,1,1)
    DECrand = np.arccos(cosdec)*180./np.pi-90.

    pix=hp.ang2pix(nside,RArand,DECrand,lonlat=True)

    if(qso_ang_mask[pix]==1):

        RAran[nr]=RArand
        DECran[nr]=DECrand

        nr=nr+1

hp.mollview(qso_ang_mask,title='QSO angular mask and random angular positions')
hp.projscatter(RAran,DECran,lonlat=True,s=0.1)
plt.show()
#---------------------------------------------------------

# RANDOM REDSHIFT DISTRIBUTION FROM THE QSO REDSHIFT DISTRIBUTION____________________________

# function to find y knowing x............................
def find_x(x,y,y0,nbins):
    for i in range(1,nbins-1):
        x0=-99.
        if((y[i-1]<y0)&(y[i+1]>y0)):
            x0 = x[i]
            break
    return(x0)
#---------------------------------------------------------

# CUMULATIVE DISTRIBUTION FUNCTION OF THE QSO REDSHIFT DISTRIBUTION
Nbins=1000

count, bins_count = np.histogram(D['z'], bins=Nbins)
pdf = count / sum(count)
cdf = np.cumsum(pdf)

plt.plot(bins_count[0:Nbins],cdf,label='cumulative distribution function of the QSO sample')
plt.legend()
plt.xlabel('z')
plt.show()

# USE OF THE CDF TO SELECT RANDOM REDSHIFTS FOLLOWING THE QSO REDSHIFT DISTRIBUTION

Zran=np.zeros(Nran)
np.random.seed()

for i in range(Nran):

    rndY=np.random.random()

    Zran[i]=find_x(bins_count[0:Nbins],cdf,rndY,Nbins)
    while(Zran[i]==-99.):
        rndY=np.random.random()
        Zran[i]=find_x(bins_count[0:Nbins],cdf,rndY,Nbins)
        

plt.hist(D['z'],density=True,bins=25,histtype='step',label='QSO z distribution') 
plt.hist(Zran,density=True,bins=25,histtype='step',label='Random z distribution')
plt.legend()
plt.xlabel('z')
plt.show() 
#---------------------------------------------------------




