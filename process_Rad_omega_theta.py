import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
import scipy.integrate as integrate
import scipy.special as special 

######## Constant defined here ########
pi        =     3.1415926535897932384626
pi2d      =     180./pi
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
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 20,  
       }  

###############read data into array################
from_path = './C-Rad/'
grid_omega_x  = np.loadtxt(from_path+'grid_omega_x.txt')
grid_theta_y  = np.loadtxt(from_path+'grid_theta_y.txt')
grid_phi_z    = np.loadtxt(from_path+'grid_phi_z.txt')
data_I        = np.loadtxt(from_path+'data.txt')

grid_theta_y = pi-grid_theta_y

print(grid_omega_x.size)
print(grid_theta_y.size)
print(grid_phi_z.size)
print(data_I.size)

data_I  =  data_I.reshape(grid_omega_x.size, grid_theta_y.size, grid_phi_z.size)
print(data_I.shape)

data_omega_theta = np.sum(data_I,axis=2)

X, Y = np.meshgrid(grid_omega_x, grid_theta_y*pi2d)

levels = np.linspace(0, 4, 51)

#norm_x = matplotlib.colors.Normalize()

data_omega_theta = np.log10(data_omega_theta+1.)

plt.contourf(X, Y, data_omega_theta.T, levels=levels, cmap='magma')
cbar=plt.colorbar(ticks=np.linspace(0.0, 4, 5))
cbar.set_label(r'$log_{10}\frac{dI}{\sin\theta d\theta d\omega}$'+' [A.U.]', fontdict=font)
cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
#plt.plot(x_omega,np.sum(np.sum(data_I_t,axis=0),axis=0),'-b',linewidth=3)
#### manifesting colorbar, changing label and axis properties ####
#plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
plt.ylabel(r'$\theta$'+' [degree]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.xscale('log')
#plt.yscale('log')
#plt.xlim(0,400)
#plt.legend(loc='upper right',fontsize=16,framealpha=1.0)
#plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
#plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,wspace=None, hspace=None)

fig = plt.gcf()
fig.set_size_inches(10.0, 8.5)
fig.savefig('./spectral_omega_theta.png',format='png',dpi=160)
plt.close("all")

