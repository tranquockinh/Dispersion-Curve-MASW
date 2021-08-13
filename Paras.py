import numpy as np
import matplotlib.pyplot as plt

# Initializations
fs = 1000 ## Hz
N = 24 ## number of geophones
x1 = 10 ## nearest-offset
dx = 1 ## geophone spacing
du = 1/75 ## Scale factor for offset between traces 
# Forward set of seismix test
LoadData = np.loadtxt('SampleData.dat', dtype='float', delimiter=None)

DataCut = LoadData[:,:N]
u = np.zeros((np.size(DataCut,0), N), dtype=float)

for i in range(np.size(DataCut[:,0])):
    for j in range(np.size(DataCut[0,:])):        

        u[i,j] = DataCut[i,j]

Tmax = (1/fs) * (len(u[:,0]) - 1) ## total recording period [s]
T = np.linspace(0,Tmax,len(u[:,0])) ## time interval array
L = (N - 1) * dx
x = np.arange(x1, L + x1 + 1.0E-5, dx) ## start from nearest-offset x1

FigFontSize = 10

# Dispersion parameters
cT_min = 50 ## [m/s]
cT_max = 400 ## [m/s]
delta_cT = 1 ## [m/s]
