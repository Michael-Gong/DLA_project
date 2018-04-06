%matplotlib inline
#import sdf
import matplotlib
import matplotlib as mpl
mpl.style.use('https://raw.githubusercontent.com/Michael-Gong/DLA_project/master/style')
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
from mpl_toolkits.mplot3d import Axes3D
import random
from mpl_toolkits import mplot3d
from matplotlib import rc
import matplotlib.transforms as mtransforms
import sys
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
font = {'family' : 'Carlito',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 25,
       }

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

def momentum_to_theta(arg_x, arg_y):
    if arg_x >= 0:
        return np.arctan(arg_y/arg_x)
    elif arg_y > 0:
        return np.arctan(arg_y/arg_x)+np.pi
    else:
        return np.arctan(arg_y/arg_x)-np.pi





part_number=1
nsteps=20001

insert1='./Data/'
insert_n='_0'
t1=np.loadtxt(insert1+'t'+insert_n+'.txt')
z1=np.loadtxt(insert1+'z'+insert_n+'.txt')
y1=np.loadtxt(insert1+'y'+insert_n+'.txt')
x1=np.loadtxt(insert1+'x'+insert_n+'.txt')
px1=np.loadtxt(insert1+'px'+insert_n+'.txt')
py1=np.loadtxt(insert1+'py'+insert_n+'.txt')
pz1=np.loadtxt(insert1+'pz'+insert_n+'.txt')
#ey=np.loadtxt(insert+'e_part'+'.txt')
#bz=np.loadtxt(insert+'b_part'+'.txt')
#ay=np.loadtxt(insert+'a_part'+'.txt')
radn1=np.loadtxt(insert1+'radn'+insert_n+'.txt')
radt1=np.loadtxt(insert1+'radt'+insert_n+'.txt')
radpx1=np.loadtxt(insert1+'rad_px'+insert_n+'.txt')
opt1=np.loadtxt(insert1+'opt'+insert_n+'.txt')
eta1=np.loadtxt(insert1+'eta'+insert_n+'.txt')

t=np.reshape(t1,(part_number,nsteps))
x=np.reshape(x1,(part_number,nsteps))
y=np.reshape(y1,(part_number,nsteps))
z=np.reshape(z1,(part_number,nsteps))
px=np.reshape(px1,(part_number,nsteps))
py=np.reshape(py1,(part_number,nsteps))
pz=np.reshape(pz1,(part_number,nsteps))
#ey=np.reshape(ey,(part_number,nsteps))
#ay=np.reshape(ay,(part_number,nsteps))
radn=np.reshape(radn1,(part_number,nsteps))
radt=np.reshape(radt1,(part_number,nsteps))
radpx=np.reshape(radpx1,(part_number,nsteps))
opt=np.reshape(opt1,(part_number,nsteps))
eta=np.reshape(eta1,(part_number,nsteps))

gamma=np.sqrt(px**2+py**2+1)

R_dep=gamma-px

py_0 = 1.00*int(150)
R_max = py_0-(radt-radpx)

y_max0=(py_0/0.02)**0.5
y_max=((py_0-(radt-radpx))/0.02)**0.5/2.0/np.pi

index=0

radn_x=(radn[index,1:]-radn[index,:-1])
radt_x=(radt[index,1:]-radt[index,:-1])

condition = np.where(radt_x>2)


#momentum_to_theta(arg_x, arg_y)

#index2=np.where(abs(y[index,:]/2/np.pi-0.0) < 1.0e-3)
plt.subplot(5,1,1)
plt.scatter(x[index,:]/2/np.pi, y[index,:]/2/np.pi, c=R_dep[index,:], s=4, cmap='magma_r', edgecolors='None')
#plt.scatter((t[index,:]-x[index,:])/2/np.pi, y_max, c=R_dep[index,:], s=1, cmap='rainbow', edgecolors='None')
plt.plot(x[index,:]/2/np.pi, y_max[index,:],'--k',linewidth=2.5,label=r'$Reduced\ w=1$')
plt.plot(x[index,:]/2/np.pi, -y_max[index,:],'--k',linewidth=2.5)
#plt.legend(loc='upper right')
#plt.colorbar()
cbar=plt.colorbar(ticks=np.linspace(np.min(R_dep[index,:]), np.max(R_dep[index,:]), 5),shrink=1)
cbar.set_label(r'$R$', fontdict=font)

plt.scatter(x[index,condition]/2/np.pi,y[index,condition]/2/np.pi,c=(radt_x[condition])[np.newaxis,:], s=50, cmap='rainbow', edgecolors='None')
cbar=plt.colorbar(ticks=np.linspace(np.min(radt_x[radn_x>0]), np.max(radt_x[radn_x>0]), 5), shrink=1)
cbar.set_label(r'$photon\ energy\ [m_ec^2]$', fontdict=font)
#plt.xlim(47,53)
#plt.ylim(-1.1,-0.7)
plt.xlabel(r'$x\ [\mu m]$',fontdict=font)
plt.ylabel(r'$y\ [\mu m]$',fontdict=font)
plt.xticks(fontsize=30); plt.yticks(fontsize=30);
plt.title(r'$electron\ for\ a_0=200,\ \alpha=0.02,\  p_0=150\ and\ simulation\ time\ 1000T_0$',fontdict=font)

arg_px = px[index,condition]
arg_py = py[index,condition]
radt_x=(radt[index,1:]-radt[index,:-1])
arg_gg = radt_x[condition]
theta_x = np.zeros_like(arg_px)
#print(arg_py/arg_px)
for i in range(np.size(theta_x)):
    theta_x[index,i]=momentum_to_theta(arg_px[index,i],arg_py[index,i])

plt.subplot(5,1,2)
plt.scatter(theta_x/np.pi*180, arg_gg, c=np.linspace(1,np.size(theta_x),np.size(theta_x))[np.newaxis,:], s=20, cmap='nipy_spectral', edgecolors='None')
cbar=plt.colorbar(ticks=np.linspace(1, np.size(theta_x), 5), shrink=1)# orientation='horizontal', shrink=0.2)
cbar.set_label(r'$Nth$', fontdict=font)
plt.xlim(-50,50)
#print(theta_x)
plt.xlabel(r'$\theta\ [degree]$',fontdict=font)
plt.ylabel(r'$\gamma$',fontdict=font)
plt.xticks(fontsize=30); plt.yticks(fontsize=30);
#plt.ylim(0,2000.0)


theta_grid = np.linspace(-180.0, 180.0, 1001)
theta_energy = np.zeros_like(theta_grid)
theta_x/np.pi*180
for i in range(np.size(theta_x/np.pi*180)):
    theta_energy[np.min(np.where(theta_x[index,i]/np.pi*180 < theta_grid))]+=arg_gg[i]
#print(theta_energy)
plt.subplot(5,1,3)
plt.plot(theta_grid, theta_energy, '-r', linewidth=3)
plt.xlim(-20,20)
#print(theta_x)
plt.xlabel(r'$\theta\ [degree]$',fontdict=font)
plt.ylabel(r'$\gamma$',fontdict=font)
plt.xticks(fontsize=30); plt.yticks(fontsize=30);


plt.subplot(5,1,4)
plt.scatter(px[index,:],py[index,:], c=R_dep[index,:], s=4, cmap='magma_r', edgecolors='None')
#plt.plot((py[index,:]),px[index,:],'--k',linewidth=2.5,label='No RR')
#plt.legend(loc='upper right')
#plt.colorbar()
cbar=plt.colorbar(ticks=np.linspace(np.min(R_dep[index,:]), np.max(R_dep[index,:]), 5))
cbar.set_label(r'$R$', fontdict=font)
#plt.xlim(47,53)
#plt.ylim(-1.1,-0.7)
#plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)
plt.scatter(px[index,condition],py[index,condition],c=(radt_x[condition])[np.newaxis,:], s=50, cmap='rainbow', edgecolors='None')
cbar=plt.colorbar(ticks=np.linspace(np.min(radt_x[radn_x>0]), np.max(radt_x[radn_x>0]), 5), shrink=1)
cbar.set_label(r'$photon\ energy\ [m_ec^2]$', fontdict=font)
plt.xlabel('$P_x\ [m_ec]$',fontdict=font)
plt.ylabel('$P_y\ [m_ec]$',fontdict=font)
plt.xticks(fontsize=30); plt.yticks(fontsize=30);


plt.subplot(5,1,5)
#plt.scatter(x[index,condition]/2/np.pi,y[index,condition]/2/np.pi,c=(radt_x[condition])[np.newaxis,:], s=50, cmap='rainbow', edgecolors='None')
#cbar=plt.colorbar(ticks=np.linspace(np.min(radt_x[radn_x>0]), np.max(radt_x[radn_x>0]), 5), shrink=1)
#cbar.set_label(r'$photon\ energy\ [m_ec^2]$', fontdict=font)

plt.scatter(x[index,:]/2/np.pi, y[index,:]/2/np.pi, c=eta[index,:], s=4, cmap='terrain', edgecolors='None')
#plt.scatter((t[index,:]-x[index,:])/2/np.pi, y_max, c=R_dep[index,:], s=1, cmap='rainbow', edgecolors='None')
plt.plot(x[index,:]/2/np.pi, y_max[index,:],'--k',linewidth=2.5,label=r'$Reduced\ w=1$')
plt.plot(x[index,:]/2/np.pi, -y_max[index,:],'--k',linewidth=2.5)
#plt.legend(loc='upper right')
#plt.colorbar()
cbar=plt.colorbar(ticks=np.linspace(np.min(eta[index,:]), np.max(eta[index,:]), 5),shrink=1)
cbar.set_label(r'$\chi_e$', fontdict=font)

#plt.xlim(47,53)
#plt.ylim(-1.1,-0.7)
plt.xlabel(r'$x\ [\mu m]$',fontdict=font)
plt.ylabel(r'$y\ [\mu m]$',fontdict=font)
plt.xticks(fontsize=30); plt.yticks(fontsize=30);
#plt.title(r'$electron\ for\ a_0=200,\ \alpha=0.02,\  p_0=150$',fontdict=font)

plt.subplots_adjust(top=0.95, bottom=0.10, left=0.15, right=0.95, hspace=0.25, wspace=0.30)

fig = plt.gcf()
#fig.set_size_inches(30, 15)
fig.set_size_inches(15, 33.5)
fig.savefig('./qed_photon.png',format='png',dpi=160)
#plt.close("all")
