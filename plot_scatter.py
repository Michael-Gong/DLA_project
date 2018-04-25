#%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
font = {'family' : 'monospace',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 20,
        }
part_number=100
nsteps=20002
insert='./Data/'
t=np.loadtxt(insert+'t'+'.txt')
y=np.loadtxt(insert+'y'+'.txt')
x=np.loadtxt(insert+'x'+'.txt')
px=np.loadtxt(insert+'px'+'.txt')
py=np.loadtxt(insert+'py'+'.txt')
#ey=np.loadtxt(insert+'e_part'+'.txt')
#bz=np.loadtxt(insert+'b_part'+'.txt')
#ay=np.loadtxt(insert+'a_part'+'.txt')
radn=np.loadtxt(insert+'radn'+'.txt')
radt=np.loadtxt(insert+'radt'+'.txt')
opt=np.loadtxt(insert+'opt'+'.txt')
eta=np.loadtxt(insert+'eta'+'.txt')

t=np.reshape(t,(part_number,nsteps))
x=np.reshape(x,(part_number,nsteps))
y=np.reshape(y,(part_number,nsteps))
px=np.reshape(px,(part_number,nsteps))
py=np.reshape(py,(part_number,nsteps))
#ey=np.reshape(ey,(part_number,nsteps))
#ay=np.reshape(ay,(part_number,nsteps))
radn=np.reshape(radn,(part_number,nsteps))
radt=np.reshape(radt,(part_number,nsteps))
opt=np.reshape(opt,(part_number,nsteps))
eta=np.reshape(eta,(part_number,nsteps))

print(np.where(py[:,-1] > 0))

gamma=np.sqrt(px**2+py**2+1)
    
rrn=(radn[:,1:-1]-radn[:,0:-2])
rrt=(radt[:,1:-1]-radt[:,0:-2])

series=np.where(py[:,-1] > 10)

for index in np.reshape(series,(np.size(series),)):
    which=np.where(rrt[index,:] >= 1)
    plt.subplot(7,1,1)
    #print(rrt.shape,x[index,which].shape)
    plt.scatter(x[index,which]/2/np.pi, y[index,which]/2/np.pi, c=rrt[index,which], s=500, cmap='rainbow', edgecolors='None')
    plt.plot(x[index,:]/2/np.pi,y[index,:]/2/np.pi,'--k',linewidth=2.5,label='QED RR')
    plt.legend(loc='upper right')
    plt.colorbar()
    plt.xlim(48,52)
    plt.ylim(-1.1,-0.7)
    plt.xlabel('$x [\mu m]$',fontdict=font)
    plt.ylabel('$y [\mu m]$',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)

    plt.subplot(7,1,2)
    plt.scatter(x[index,which]/2/np.pi, px[index,which], c=rrt[index,which], s=500, cmap='rainbow', edgecolors='None')
    plt.plot(x[index,:]/2/np.pi,px[index,:],'--k',linewidth=2.5,label='QED RR')
    plt.legend(loc='upper right')
    plt.colorbar()
    plt.xlim(48,52)
    #plt.ylim(-1.1,-0.7)
    plt.xlabel('$x [\mu m]$',fontdict=font)
    plt.ylabel('$px [\mu m]$',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)
    
    plt.subplot(7,1,3)
    plt.scatter(x[index,which]/2/np.pi, py[index,which], c=rrt[index,which], s=500, cmap='rainbow', edgecolors='None')
    plt.plot(x[index,:]/2/np.pi,py[index,:],'--k',linewidth=2.5,label='QED RR')
    plt.legend(loc='upper right')
    plt.colorbar()
    plt.xlim(48,52)
    #plt.ylim(-1.1,-0.7)
    plt.xlabel('$x [\mu m]$',fontdict=font)
    plt.ylabel('$py [\mu m]$',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)
    
    plt.subplot(7,1,4)
    plt.scatter(x[index,which]/2/np.pi, radn[index,which], c=rrt[index,which], s=500, cmap='rainbow', edgecolors='None')
    plt.plot(x[index,:]/2/np.pi,radn[index,:],'-r',linewidth=2.5,label='QED RR')
    plt.legend(loc='upper right')
    plt.colorbar()
    plt.xlim(48,52)
    #plt.ylim(-1.5,-0.5)
    plt.xlabel('$x [\mu m]$',fontdict=font)
    plt.ylabel('$N [count]$',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)

    plt.subplot(7,1,5)
    plt.scatter(x[index,which]/2/np.pi, radt[index,which], c=rrt[index,which], s=500, cmap='rainbow', edgecolors='None')
    plt.plot(x[index,:]/2/np.pi,radt[index,:],'-r',linewidth=2.5,label='QED RR')
    plt.legend(loc='upper right')
    plt.colorbar()
    plt.xlim(48,52)
    #plt.ylim(-1.5,-0.5)
    plt.xlabel('$x [\mu m]$',fontdict=font)
    plt.ylabel('$\gamma$',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)

    plt.subplot(7,1,6)
    plt.scatter(x[index,which]/2/np.pi, eta[index,which], c=rrt[index,which], s=500, cmap='rainbow', edgecolors='None')
    plt.plot(x[index,:]/2/np.pi,eta[index,:],'-r',linewidth=2.5,label='QED RR')
    plt.legend(loc='upper right')
    plt.colorbar()
    plt.xlim(48,52)
    #plt.ylim(-1.5,-0.5)
    plt.xlabel('$x [\mu m]$',fontdict=font)
    plt.ylabel('$eta$',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)

    plt.subplot(7,1,7)
    plt.scatter(x[index,which]/2/np.pi, opt[index,which], c=rrt[index,which], s=500, cmap='rainbow', edgecolors='None')
    plt.plot(x[index,:]/2/np.pi,opt[index,:],'-r',linewidth=2.5,label='QED RR')
    plt.legend(loc='upper right')
    plt.colorbar()
    plt.xlim(48,52)
    #plt.ylim(-1.5,-0.5)
    plt.xlabel('$x [\mu m]$',fontdict=font)
    plt.ylabel('$opt$',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)

    #plt.show()
    #lt.figure(figsize=(100,100))
    fig = plt.gcf()
    fig.set_size_inches(10, 50)
    fig.savefig('./limit1/'+str(index).zfill(4)+'.png',format='png',dpi=60)
    plt.close("all")
    #np.max(radn[index,1:-1]-radn[index,0:-2])
    #np.max(radt[index,1:-1]-radt[index,0:-2])
    #print(rrt)
