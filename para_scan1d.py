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
nsteps=200007
axis_b=np.linspace(0,100,101)
#axis_w=np.linspace(30,90,16)
result1=np.zeros([101])
#rebo2=np.zeros([16,16])
insert1='./'
#insert2='epoch2dqe/'
insert_n='_0'
a0=1.0
for ib in range(101):
  px1=np.loadtxt(insert1+'Datab'+str(int(axis_b[ib]))+'/px'+insert_n+'.txt')
  py1=np.loadtxt(insert1+'Datab'+str(int(axis_b[ib]))+'/py'+insert_n+'.txt')
  px=np.reshape(px1,(part_number,nsteps))
  py=np.reshape(py1,(part_number,nsteps))
  gamma=np.sqrt(px**2+py**2+1)
  index=0
  result1[ib]=max(gamma[index,:])/(a0**2/2.0+1)

plt.plot(axis_b/100.0,result1,'-r',linewidth=2.5,label=r'$\frac{\gamma_{max}}{\gamma_{b_0=0}}$')
plt.legend(loc='upper left')
#plt.ylim(-1.1,-0.7)
plt.xlabel(r'$b_0$',fontdict=font)
plt.ylabel(r'$\frac{\gamma_{max}}{\gamma_{b_0=0}}$',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)


fig = plt.gcf()
fig.set_size_inches(8, 6)
#fig.set_size_inches(5, 4.5)
fig.savefig('1dscan_for_b.png',format='png',dpi=160)
plt.close("all")
