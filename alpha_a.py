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

#plt.scatter(theta_x/np.pi*180, arg_gg, c=np.linspace(1,np.size(theta_x),np.size(theta_x))[np.newaxis,:], s=20, cmap='nipy_spectral', edgecolors='None')
#cbar=plt.colorbar(ticks=np.linspace(1, np.size(theta_x), 5), shrink=1)# orientation='horizontal', shrink=0.2)
#cbar.set_label(r'$Nth$', fontdict=font)
#plt.xlim(-45,45)
##print(theta_x)
#plt.xlabel(r'$\theta\ [degree]$',fontdict=font)
#plt.ylabel(r'$\gamma$',fontdict=font)
##plt.xticks(fontsize=30); plt.yticks(fontsize=30);
##plt.ylim(0,2000.0)

a0=np.linspace(10,210,1001)
#alpha=0.04**1.5*a0/(4.6**0.75)
alpha= (179.0**0.5*a0**2/2.3e6-9.6*a0**2/2.03e6-1.3e1/2.03e6)**0.5
#plt.plot(a0,alpha,'-k',linewidth=4)
plt.plot(a0,(a0**2-6.5)**0.5/1000.0,'-k',linewidth=4)
alpha=0.04**1.5*a0/(4.6**0.75)
#plt.plot(a0,alpha,'--b',linewidth=4)


u =  1.0/12.5
a0_1=np.array([10,25,50,75,100,125,150,200])
alpha_1=np.array([-2+2*u,-2+6*u,-2+10*u,-2+11*u,-1+1.5*u,-1+3*u,-1+4*u,-1+5*u])
plt.scatter(a0_1,10**(alpha_1-0.25*u),marker='+',s=40,color='r')

plt.xlabel(r'$a_0$',fontdict=font)
plt.ylabel(r'$\alpha$',fontdict=font)
plt.xticks(fontsize=30); plt.yticks(fontsize=30);
plt.yscale('log')
plt.ylim(10**-2,10**0)

fig = plt.gcf()
#fig.set_size_inches(30, 15)
fig.set_size_inches(8, 4)
#fig.savefig('./bunch_theta_en.png',format='png',dpi=160)
#plt.close("all")
