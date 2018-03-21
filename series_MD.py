from scipy.integrate import odeint
#%matplotlib inline
#import sdf
import matplotlib
import matplotlib as mpl
#mpl.style.use('https://raw.githubusercontent.com/Michael-Gong/DLA_project/master/style')
matplotlib.use('agg')
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
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

font = {'family' : 'helvetica',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 25,
        }

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# function that returns dz/dt
def model(z,t,a0):
    u = z[0]
    w = z[1]
    dudt = (alpha**0.5/(C1**1.5))/(1-w**2)*a0+w*u**2/(1-w**2)-(alpha/(C1**3))*w/(1-w**2)**3+(alpha/C1)*w/(1-w**2)#(-x + u)/2.0
    dwdt = u#(-y + x)/5.0
    dzdt = [dudt,dwdt]
    return dzdt

# initial condition
P0=100.0
alpha=0.01
C1=(P0**2+1.0)**0.5
#z0 = [0,0]
z0 = [(alpha/C1)**0.5,0]

# number of time points
nsteps=300000

# time points
t = np.linspace(0,30*2*np.pi,nsteps)

# step input laser a0
a0 = -50.0*np.cos(t)

# store solution
u = np.empty_like(t)
w = np.empty_like(t)
# record initial conditions
u[0] = z0[0]
w[0] = z0[1]

# solve ODE
for i in range(1,nsteps):
    # span for next time step
    tspan = [t[i-1],t[i]]
    # solve for next step
    z = odeint(model,z0,tspan,args=(a0[i],))
    # store solution for plotting
    u[i] = z[1][0]
    w[i] = z[1][1]
    # next initial condition
    z0 = z[1] 

R=C1-alpha*w**2*C1/alpha    
py=u*R*(C1/alpha)**0.5

gamma=(1+py**2+R**2)/2/R
px=gamma-R
q=w*(C1/alpha)**0.5

d_work=py/R*a0

term1=alpha*q*(py/R)**2/R
term2=-alpha*q/R**3
term3=alpha*q/R
term4=a0/R

n_min=0 # 104999 #99999# 74999 #35000
n_max=299999 # 124999 #5000 #399999
lgR=R[n_min:n_max]
#lgR=np.log10(R[n_min:n_max])
# plot results

x_grid = np.linspace(-1.5,31.5,601)
y_grid = np.linspace(-1.2,1.2,201)

[x_grid,y_grid] = np.meshgrid(x_grid,y_grid)
Ey = np.cos(2*np.pi*x_grid)

eee=np.max([-np.min(Ey),np.max(Ey)])
levels = np.linspace(-eee, eee, 64)

length_x = 2.0

for n_snap in range(0, 50000, 500):
        print('This is '+str(n_snap).zfill(5)+' of '+str(20000))
        start_x = -1.0+(t[-1]-t[-2])*n_snap/2/np.pi
        ax=plt.subplot(2,2,1)
        #### manifesting colorbar, changing label and axis properties ####
        #cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee])
        #cbar.set_label('Normalized electric field',fontdict=font)        
        plt.contourf(x_grid, y_grid, Ey, levels=levels, cmap=cm.bwr, alpha=.7)
        plt.scatter(t[n_min:n_max]/2/np.pi,w[n_min:n_max], c=lgR, s=5, cmap='rainbow', edgecolors='None')
        plt.scatter(t[n_min+n_snap]/2/np.pi,w[n_min+n_snap], color='r', marker='o', s=200, edgecolors='k')
        #plt.plot((t[index,:])/2/np.pi,np.sqrt(px[index,:]**2+py[index,:]**2+1),'--k',linewidth=2.5,label='No RR')
        #plt.legend(loc='upper right')
        #cbar=plt.colorbar(ticks=np.linspace(np.min(gamma), np.max(gamma), 5))
        #cbar=plt.colorbar(ticks=np.linspace(np.min(lgR), np.max(lgR), 5))
        #cbar.set_label(r'$R$', fontdict=font)#plt.xlim(47,53)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        #plt.xlabel(r'$\xi\ [2\pi]$',fontdict=font)
        plt.ylabel(r'$w\ [\sqrt{\frac{\alpha}{C_1}}\frac{\lambda}{2\pi}]$',fontdict=font)
        plt.xticks(fontsize=1.0); plt.yticks(fontsize=26);
        plt.ylim(-1.025,1.025)
        plt.xlim(start_x,start_x+length_x)
        #plt.legend(loc='best')


        ax=plt.subplot(2,2,2)
        #plt.plot(t,a0,'g:',label=r'$a_0(t)$')
        #plt.plot(t,u,'b-',label='u(t)')
        #plt.plot(t/2/np.pi,w,'r-',label='w(t)')
        #n_min=0 # 104999 #99999# 74999 #35000
        #n_max=299999 # 124999 #5000 #399999
        x1=t[n_min:n_max]/2/np.pi
        y1=d_work[n_min:n_max]
        y2=np.zeros_like(x1)
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        #ax.fill_between(x1, 0, 1, where=abs(x1-11.5)<0.25 , facecolor='cyan', alpha=0.2, transform=trans)
        ax.fill_between(x1,y1,y2,where=y1>= y2, facecolor='green', alpha=0.3, interpolate=True)
        ax.fill_between(x1,y1,y2,where=y1<= y2, facecolor='red', alpha=0.3, interpolate=True)
        #ax.fill(t[n_min:n_max]/2/np.pi,d_work[n_min:n_max], 'b', alpha=0.3)
        #plt.scatter(t[n_min:n_max]/2/np.pi,d_work[n_min:n_max], c=lgR, s=1, cmap='rainbow', edgecolors='None')
        plt.scatter(t[n_min:n_max]/2/np.pi,d_work[n_min:n_max], c=lgR, s=5, cmap='rainbow', edgecolors='None')
        plt.scatter(t[n_min+n_snap]/2/np.pi,d_work[n_min+n_snap], color='r', marker='o', s=200, edgecolors='k')
        #plt.plot(t[n_min:n_max]/2/np.pi,d_work[n_min:n_max],'-k',linewidth=2)
        #plt.plot((t[index,:])/2/np.pi,np.sqrt(px[index,:]**2+py[index,:]**2+1),'--k',linewidth=2.5,label='No RR')
        #plt.legend(loc='upper right')
        #cbar=plt.colorbar(ticks=np.linspace(np.min(gamma), np.max(gamma), 5))
        #cbar=plt.colorbar(ticks=np.linspace(np.min(lgR), np.max(lgR), 5))
        #cbar.set_label(r'$log_{10}R$', fontdict=font)#plt.xlim(47,53)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        #plt.xlabel(r'$\xi\ [2\pi]$',fontdict=font)
        plt.ylabel(r'$\frac{dE_k}{d\xi}$',fontdict=font)
        plt.grid(color='k', linestyle='--', linewidth=1)
        plt.xticks(fontsize=1); plt.yticks(fontsize=26);
        plt.ylim(-709.9,709.9)
        #plt.xlim(11.0,12.0)
        plt.xlim(start_x,start_x+length_x)
        #plt.legend(loc='best')

        ax=plt.subplot(2,2,3)
        #plt.plot(t,a0,'g:',label=r'$a_0(t)$')
        #plt.plot(t,u,'b-',label='u(t)')
        #plt.plot(t/2/np.pi,gamma,'k-',label=r'$\gamma$')
        plt.scatter(t[n_min:n_max]/2/np.pi,gamma[n_min:n_max], c=lgR, s=5, cmap='rainbow', edgecolors='None')
        plt.scatter(t[n_min+n_snap]/2/np.pi,gamma[n_min+n_snap], color='r', marker='o', s=200, edgecolors='k')
        #plt.plot((t[index,:])/2/np.pi,np.sqrt(px[index,:]**2+py[index,:]**2+1),'--k',linewidth=2.5,label='No RR')
        #plt.legend(loc='upper right')
        #cbar=plt.colorbar(ticks=np.linspace(np.min(gamma), np.max(gamma), 5))
        #cbar.set_label(r'$\gamma$', fontdict=font)#plt.xlim(47,53)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        plt.xlabel(r'$\xi\ [2\pi]$',fontdict=font)
        plt.ylabel(r'$\gamma$',fontdict=font)
        plt.grid(color='k', linestyle='--', linewidth=1)
        plt.xticks(fontsize=26); plt.yticks(fontsize=26);
        #plt.ylim(-1.1,1.1)
        plt.ylim(-109.0, 2409.0)
        #plt.xlim(11.0,12.0)
        plt.xlim(start_x,start_x+length_x)
        #plt.legend(loc='best')

        
        ax=plt.subplot(2,2,4)
        #plt.plot(t,a0,'g:',label=r'$a_0(t)$')
        #plt.plot(t,u,'b-',label='u(t)')
        #plt.plot(px,py,'g-',label=r'$p_x-p_y$')
        p_py = np.linspace(-750.0,750.0,401)
        p_px = np.linspace(-50.0,2650.0,601)

        [p_px,p_py] = np.meshgrid(p_px,p_py)

        p_R = (1+p_px**2+p_py**2)**0.5-p_px
        levels = np.linspace(np.min(p_R), np.max(p_R), 32)
        plt.contour(p_px, p_py, p_R, np.array([100.0]),linestyles='dashed',linewidth=0.5, alpha=0.35)
        #plt.scatter(px2[n_min2:n_max2],py2[n_min2:n_max2], c=lgR2, s=3, cmap='plasma_r', edgecolors='None')
        plt.scatter(px[n_min:n_max],py[n_min:n_max], c=lgR, s=3, cmap='rainbow', edgecolors='None')
        plt.scatter(px[n_min+n_snap],py[n_min+n_snap], color='r', marker='o', s=200, edgecolors='k')
        #cbar=plt.colorbar(ticks=np.linspace(np.min(lgR), np.max(lgR), 5))
        #cbar.set_label(r'$R$', fontdict=font)#plt.xlim(47,53)
        #plt.scatter(px1[n_min:n_max],py1[n_min:n_max], s=3, color='k', edgecolors='None')
        #plt.plot((t[index,:])/2/np.pi,np.sqrt(px[index,:]**2+py[index,:]**2+1),'--k',linewidth=2.5,label='No RR')
        #plt.legend(loc='upper right')
        #cbar=plt.colorbar(ticks=np.linspace(np.min(gamma), np.max(gamma), 5))

        #cbar.set_label(r'$\gamma$', fontdict=font)#plt.xlim(47,53)
        plt.xlabel(r'$p_x\ [m_ec]$',fontdict=font)
        plt.ylabel(r'$p_y\ [m_ec]$',fontdict=font)
        plt.xticks(fontsize=30); plt.yticks(fontsize=30);
        plt.ylim(-745.0,745.0)
        plt.xlim(-99.0,2999.0)
        #plt.legend(loc='best')

        #plt.show()
        plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.05, wspace=0.35)
        fig = plt.gcf()
        fig.set_size_inches(18, 14)
        #fig.set_size_inches(5, 4.5)
        fig.savefig('./series_MD/series'+str(n_snap).zfill(5)+'.png',format='png',dpi=160)
        plt.close("all")
