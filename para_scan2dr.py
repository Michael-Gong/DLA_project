#%matplotlib inline
#import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
#from colour import Color

######## Constant defined here ########
pi        =     3.1415926535897932384626
q0        =     1.602176565e-19 # C
m0        =     9.10938291e-31  # kg
v0        =     2.99792458e8    # m/s^2
kb        =     1.3806488e-23   # J/K
mu0       =     4.0e-7*pi       # N/A^2
epsilon0  =     8.8541878176203899e-12 # F/m
h_planck  =     6.62606957e-34  # J s
wavelength=     1.0e-6
frequency =     v0*2*pi/wavelength

exunit    =     m0*v0*frequency/q0
bxunit    =     m0*frequency/q0
denunit    =     frequency**2*epsilon0*m0/q0**2
#print 'electric field unit: '+str(exunit)
#print 'magnetic field unit: '+str(bxunit)
#print 'density unit nc: '+str(denunit)

font = {'family' : 'monospace',  
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 20,  
       }  



part_number=1
nsteps=100004
axis_r=np.linspace(0,160,41)
axis_a=np.linspace(0,800,41)
#axis_w=np.linspace(30,90,16)
result1=np.zeros([41,41])
#rebo2=np.zeros([16,16])
insert1='./'
#insert2='epoch2dqe/'
insert_n='_0'
for ir in range(41):
    for ia in range(41):
        px1=np.loadtxt(insert1+'Datar'+str(int(axis_r[ir]))+'a'+str(int(axis_a[ia]))+'/px'+insert_n+'.txt')
        py1=np.loadtxt(insert1+'Datar'+str(int(axis_r[ir]))+'a'+str(int(axis_a[ia]))+'/py'+insert_n+'.txt')
        px=np.reshape(px1,(part_number,nsteps))
        py=np.reshape(py1,(part_number,nsteps))
        gamma=np.sqrt(px**2+py**2+1)
        index=0
        a0=axis_a[ia]/100.0
        result1[ir,ia]=max(gamma[index,:])/(a0**2/2.0+1)
        print("finished "+str((ia*1.0+ir*41.0)/(41.0*41.0)*100)+" %")
np.savetxt('./txt/result_2d.txt', result1)
np.savetxt('./txt/axis_r.txt', axis_r/100)
np.savetxt('./txt/axis_a.txt', axis_a/100)
x=axis_r/100
y=axis_a/100
z_min=np.min(result1)
z_max=np.max(result1);
plt.imshow(result1, cmap='terrain', vmin=z_min, vmax=z_max,interpolation='none')
#plt.axis([x.min(), x.max(), y.min(), y.max()])
plt.colorbar()
plt.title(r'$\frac{\gamma_{max}}{\gamma_{b_0=0}}$')
#plt.ylim(-1.1,-0.7)
plt.xlabel(r'$b_0/a_0$',fontdict=font)
plt.ylabel(r'$a_0$',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)


fig = plt.gcf()
fig.set_size_inches(8, 6)
#fig.set_size_inches(5, 4.5)
fig.savefig('./txt/2dscan_for_r_a.png',format='png',dpi=160)
plt.close("all")
