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
axis_b=np.linspace(0,40000,51)
#axis_a=np.array([2000,5000,10000,15000])
#axis_a=np.linspace(100,2500,41)
axis_p=np.linspace(0,200,51)
#axis_w=np.linspace(30,90,16)
#enhancement=np.zeros([51,51])
#ymax=np.zeros([51,51])
#bchf=np.zeros([51,51])
R=np.zeros([51,51])

#rebo2=np.zeros([16,16])
insert1='./'
#insert2='epoch2dqe/'
insert_n='_0'
for ib in range(51):
    for ip in range(51):
        px1=np.loadtxt(insert1+'Datab'+str(int(axis_b[ib]))+'p'+str(int(axis_p[ip]))+'/px'+insert_n+'.txt')
        py1=np.loadtxt(insert1+'Datab'+str(int(axis_b[ib]))+'p'+str(int(axis_p[ip]))+'/py'+insert_n+'.txt')
       # y1=np.loadtxt(insert1+'Datab'+str(int(axis_b[ib]))+'p'+str(int(axis_p[ip]))+'/y'+insert_n+'.txt')
       # efield1=np.loadtxt(insert1+'Datab'+str(int(axis_b[ib]))+'p'+str(int(axis_p[ip]))+'/e_part'+insert_n+'.txt')
       # bfield1=np.loadtxt(insert1+'Datab'+str(int(axis_b[ib]))+'p'+str(int(axis_p[ip]))+'/b_part'+insert_n+'.txt')
        px=np.reshape(px1,(part_number,nsteps))
        py=np.reshape(py1,(part_number,nsteps))
       # grid_y=np.reshape(y1,(part_number,nsteps))
       # efield=np.reshape(efield1,(part_number,nsteps))
       # bfield=np.reshape(bfield1,(part_number,nsteps))
        gamma=np.sqrt(px**2+py**2+1)
        index=0
        #a0=axis_a[ia]/100.0
        a0=200.0
        #enhancement[ib,ip]=max(gamma[index,:])/(a0**2/2.0+1)
        #ymax[ib,ip]=max(abs(grid_y[index,:]))/2/np.pi
        #bchf[ib,ip]=max(abs(bfield[index,:]-efield[index,:]))
        R[ib,ip]=np.min(gamma[index,:]-px[index,:])
        print("finished "+str((ip*1.0+ib*51.0)/(51.0*51.0)*100)+" %")
#np.savetxt('./txt/enhance.txt', enhancement)
#np.savetxt('./txt/ymax.txt', ymax)
np.savetxt('./txt/R.txt', R)
np.savetxt('./txt/axis_b.txt', axis_b/10000.0-3.5)
np.savetxt('./txt/axis_p.txt', axis_p)
#x=axis_b/100
#y=axis_a/100
#z_min=np.min(result1)
#z_max=np.max(result1);
#plt.imshow(result1, cmap='terrain', vmin=z_min, vmax=z_max,interpolation='nearest', origin='lower')
#plt.axis([x.min(), x.max(), y.min(), y.max()])
#plt.colorbar():wq




#plt.title(r'$\frac{\gamma_{max}}{\gamma_{b_0=0}}$')
#plt.ylim(-1.1,-0.7)
#plt.xlabel(r'$\alpha$',fontdict=font)
#plt.ylabel(r'$a_0$',fontdict=font)
#plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)


#fig = plt.gcf()
#fig.set_size_inches(8, 6)
#fig.set_size_inches(5, 4.5)
#fig.savefig('./txt/2dscan_for_b_a.png',format='png',dpi=160)
#plt.close("all")
