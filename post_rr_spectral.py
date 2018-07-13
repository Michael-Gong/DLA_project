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
from_path = './Data/'
px      = np.loadtxt(from_path+'px_0.txt')
py      = np.loadtxt(from_path+'py_0.txt')
pz 		= np.loadtxt(from_path+'pz_0.txt')
gg 		= (px**2+py**2+pz**2+1.0)**0.5
grid_x  = np.loadtxt(from_path+'x_0.txt')
grid_y  = np.loadtxt(from_path+'y_0.txt')
grid_z  = np.loadtxt(from_path+'z_0.txt')
time_t  = np.loadtxt(from_path+'t_0.txt')
#part_ex = np.loadtxt(from_path+'ex_part_0')
part_ex = np.zeros_like(time_t)
part_ey = np.loadtxt(from_path+'e_part_0.txt')
part_ez = np.zeros_like(time_t)
part_bx = np.zeros_like(time_t)
part_by = np.zeros_like(time_t)
part_bz = np.loadtxt(from_path+'b_part_0.txt')

e_vx    = px/gg
e_vy    = py/gg
e_vz    = pz/gg
e_ax    = (-part_ex-(e_vy*part_bz-e_vz*part_by)-e_vx*(-e_vx*part_ex-e_vy*part_ey-e_vz*part_ez))/gg
e_ay    = (-part_ey-(e_vz*part_bx-e_vx*part_bz)-e_vy*(-e_vx*part_ex-e_vy*part_ey-e_vz*part_ez))/gg
e_az    = (-part_ez-(e_vx*part_by-e_vy*part_bx)-e_vz*(-e_vx*part_ex-e_vy*part_ey-e_vz*part_ez))/gg

##################################################

#n_dirc  = np.array([0.,0.,1.])

x_omega = np.logspace(0,6,100)
y_theta = np.linspace(1,179,90)
z_phi   = np.linspace(-179,179,180)
data_I  = np.zeros([180,90,100])
data_I_t= data_I

norm_fac=1.0/4/3.14/epsilon0*q0**2/4/3.14**2/v0*frequency/(m0*v0**2)
############ To calculate integral of the dI/domega##########
for i_phi in range(np.size(z_phi)):
  for i_theta in range(np.size(y_theta)):
    for i_omega in range(np.size(x_omega)):     
      n_dirc = np.array([np.cos(y_theta[i_theta]),np.sin(y_theta[i_theta])*np.cos(z_phi[i_phi]),np.sin(y_theta[i_theta])*np.sin(z_phi[i_phi])]) 
      dt = time_t[-1]-time_t[-2]; data_1=0; data_2=0
      for i_time in range(np.size(time_t)):
        e_vel_t = np.array([e_vx[i_time], e_vy[i_time], e_vz[i_time]])
        e_acc_t = np.array([e_ax[i_time], e_ay[i_time], e_ax[i_time]])
        e_pos_t = np.array([grid_x[i_time], grid_y[i_time], grid_z[i_time]])
        e_tim_t = time_t[i_time]
        amplitude1 = np.cross(n_dirc,np.cross(n_dirc-e_vel_t,e_acc_t))/np.dot(1-np.dot(e_vel_t,n_dirc),1-np.dot(e_vel_t,n_dirc))
        phase2     = (e_tim_t-np.dot(n_dirc,e_pos_t))*x_omega[i_omega]
        data_1     = data_1+dt*amplitude1*np.cos(phase2)
        data_2     = data_2+dt*amplitude1*np.sin(phase2)
      data_I[i_phi][i_theta][i_omega]  = np.dot(data_1,data_1)+np.dot(data_2,data_2)
      data_I_t[i_phi][i_theta][i_omega]= data_I[i_phi][i_theta][i_omega]*np.sin(y_theta[i_theta]) 
      print('i_phi:',z_phi[i_phi],'i_theta:',y_theta[i_theta],'i_omega:',x_omega[i_omega],' ; dI/dw:',data_I[i_phi][i_theta][i_omega]*norm_fac)
#############################################################
data_I = norm_fac*data_I
data_I_t = norm_fac*data_I*(y_theta[-1]-y_theta[-2])*pi2d*(z_phi[-1]-z_phi[-2])*pi2d

np.savetxt('./Data/spectrum_z_phi.txt',z_phi)
np.savetxt('./Data/spectrum_y_theta.txt',y_theta)
np.savetxt('./Data/spectrum_x_omega.txt',x_omega)
np.savetxt('./Data/spectrum_dI_dw_do.txt',data_I)
np.savetxt('./Data/spectrum_dI_dw.txt',np.sum(np.sum(data_I_t,axis=0),axis=0))

plt.plot(x_omega,np.sum(np.sum(data_I_t,axis=0),axis=0),'-b',linewidth=3)
#### manifesting colorbar, changing label and axis properties ####
plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
plt.ylabel('dI/d$\omega$ [$m_ec^2/\omega_0$]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
plt.xscale('log')
plt.yscale('log')
#plt.xlim(0,400)
#plt.legend(loc='upper right',fontsize=16,framealpha=1.0)
#plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
#plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,wspace=None, hspace=None)

fig = plt.gcf()
fig.set_size_inches(10.0, 8.5)
fig.savefig('./'+'spectral_dI_dw.png',format='png',dpi=160)
plt.close("all")

