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
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors
import scipy.ndimage as ndimage



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
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 25,  
       }  

font2 = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 25,  
       }

font_size =25

to_path   = './two/'
grid_data = np.load(to_path+'ey_vph.npy')
grid_time = np.load(to_path+'grid_time.npy')/1e-15
grid_x    = np.load(to_path+'grid_x.npy')

xx  =  np.load(to_path+'xx2d.npy')
yy  =  np.load(to_path+'yy2d.npy')
px  =  np.load(to_path+'px2d.npy')
py  =  np.load(to_path+'py2d.npy')
work_x = np.load(to_path+'workx2d.npy')
work_y = np.load(to_path+'worky2d.npy')
gg  =  (1+px**2+py**2)**0.5
  
ey  =  np.load(to_path+'ey2d.npy') 
ex  =  np.load(to_path+'ex2d.npy') 
bz  =  np.load(to_path+'bz2d.npy') 

for n in range(1):
    n = 7
    plt.subplot(4,1,1)
    plt.plot(grid_time, work_y[n,:], color='tomato', linewidth=3, linestyle='-')
    plt.plot(grid_time, work_x[n,:], color='royalblue', linewidth=3, linestyle='-')
    plt.plot(grid_time, gg[n,:], color='k', linewidth=3, linestyle='--')
#    plt.scatter(grid_time, work_y[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=30, cmap='magma', edgecolors='None', alpha=1,zorder=3)
    #  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    #  cbar.set_label('$\gamma$',fontdict=font)
    #plt.ylabel('time [fs]',fontdict=font)
    plt.ylabel('$W_{x,y}$ [$m_ec^2$]',fontdict=font)
    plt.xlim(100,290)
    plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
    plt.grid()
    #plt.xscale('log')
    #plt.xlim(-200,14000)
    #plt.ylim(-800,800)
    
    
    plt.subplot(4,1,2)
    plt.plot(grid_time,-py[n,:]/gg[n,:]*ey[n,:], color='k', linewidth=3, linestyle='-')
    plt.scatter(grid_time,-py[n,:]/gg[n,:]*ey[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1400), s=30, cmap='magma', edgecolors='None', alpha=1,zorder=3)
    #  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    #  cbar.set_label('$\gamma$',fontdict=font)
    #plt.ylabel('time [fs]',fontdict=font)
    plt.ylabel('$d\gamma/dt$',fontdict=font)
    plt.xlim(100,290)
    plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
    plt.grid()
    #plt.xscale('log')
    #plt.xlim(-200,14000)
    #plt.ylim(-800,800)
    
    
    plt.subplot(4,1,3)
    plt.plot(grid_time,yy[n,:], color='k', linewidth=3, linestyle='-')
    plt.scatter(grid_time,yy[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1400), s=30, cmap='magma', edgecolors='None', alpha=1,zorder=3)
    #  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    #  cbar.set_label('$\gamma$',fontdict=font)
    #plt.ylabel('time [fs]',fontdict=font)
    plt.ylabel('$y$ [$\mu m$]',fontdict=font)
    plt.xlim(100,290)
    plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
    plt.grid()
    #plt.xscale('log')
    #plt.xlim(-200,14000)
    #plt.ylim(-800,800)
    
    plt.subplot(4,1,4) 
    #axin1 = inset_axes(ax, width='15%', height='5%', loc='upper left')
    #axin2 = inset_axes(ax, width='15%', height='5%', loc='upper center')
    X, Y = np.meshgrid(grid_x, grid_time) 
    levels = np.linspace(-60.1, 60.1, 50)
    grid_data[grid_data < -60]=-60
    grid_data[grid_data >  60]= 60
    plt.contourf(Y.T, (30+X-v0*Y*1e-15/1e-6).T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.bwr)
    #cbar=plt.colorbar(ticks=np.linspace(-60,60,3))#,orientation="horizontal")
    #cbar.set_label('$E_y$ [$m_ec\omega/e$]', fontdict=font2)
    #for n in range(7):
    #    plt.plot(xx[n,:]+30-v0*grid_time*1e-15/1e-6, grid_time,'-k')
    #plt.contourf(X, Y, grid_data.T, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.bwr)
    #n = 3
    plt.scatter(grid_time,xx[n,:]+30-v0*grid_time*1e-15/1e-6, c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1400), s=30, cmap='magma', edgecolors='None', alpha=1,zorder=3)
    #plt.xlabel('Energy [MeV]',fontdict=font)
    plt.xlabel('time [fs]',fontdict=font)
    plt.ylabel('$\Delta=x-ct+30$ [$\mu$m]',fontdict=font)
    plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
    #plt.yscale('log')
    plt.xlim(100,290)
    plt.ylim(-5,0)
    #plt.legend(loc='best',fontsize=20,framealpha=0.5)
    plt.grid()
    
    
    plt.subplots_adjust(left=0.10, bottom=0.08, right=0.95, top=0.92,
                    wspace=0.05, hspace=0.15)
    
    fig = plt.gcf()
    fig.set_size_inches(20, 15)
    fig.savefig(to_path+'wrap_trace_'+str(n).zfill(4)+'.png',format='png',dpi=160)
    plt.close("all")
    print('ok')

