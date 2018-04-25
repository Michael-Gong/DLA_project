#%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
font = {'family' : 'monospace',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 20,
        }

part_number=10000
nsteps=20002

#for insert in ['50','100','150','200','250','300']: 
for insert in ['20','40','60','80','100','120']: 
    #t=np.loadtxt(insert+'t'+'.txt')
    y=np.loadtxt('Dataa'+insert+'/y'+'.txt')
    x=np.loadtxt('Dataa'+insert+'/x'+'.txt')
    px=np.loadtxt('Dataa'+insert+'/px'+'.txt')
    py=np.loadtxt('Dataa'+insert+'/py'+'.txt')
    #ey=np.loadtxt(insert+'e_part'+'.txt')
    #bz=np.loadtxt(insert+'b_part'+'.txt')
    #ay=np.loadtxt(insert+'a_part'+'.txt')
    radn=np.loadtxt('Dataa'+insert+'/radn'+'.txt')
    radt=np.loadtxt('Dataa'+insert+'/radt'+'.txt')
    #opt=np.loadtxt(insert+'opt'+'.txt')
    #eta=np.loadtxt(insert+'eta'+'.txt')
    
    #t=np.reshape(t,(part_number,nsteps))
    x=np.reshape(x,(part_number,nsteps))
    y=np.reshape(y,(part_number,nsteps))
    px=np.reshape(px,(part_number,nsteps))
    py=np.reshape(py,(part_number,nsteps))
    #ey=np.reshape(ey,(part_number,nsteps))
    #ay=np.reshape(ay,(part_number,nsteps))
    radn=np.reshape(radn,(part_number,nsteps))
    radt=np.reshape(radt,(part_number,nsteps))
    #opt=np.reshape(opt,(part_number,nsteps))
    #eta=np.reshape(eta,(part_number,nsteps))
    
    print(np.where(py[:,-1] > 0))
    
    gamma=np.sqrt(px**2+py**2+1)
        
    rrn=(radn[:,1:-1]-radn[:,0:-2])
    rrt=(radt[:,1:-1]-radt[:,0:-2])
    
    series=np.where(py[:,-1] > 10)
    
    
    #index=75
    rrn=(radn[:,1:-1]-radn[:,0:-2])
    rrt=(radt[:,1:-1]-radt[:,0:-2])
    #which=np.where(rrt[index,:] >= 1)
    gamma=np.sqrt(px**2+py**2+1)
    #photon_gamma=0;photon_px=0;photon_py=0;
    photon_x=np.zeros(100);photon_y=np.zeros(100)
    for ith in np.arange(100):
        photon_gamma=0;photon_px=0;photon_py=0
        index=ith
        which=np.where(rrn[index,:] >= 1)
        for i in np.reshape(which,(np.size(which),)):
            photon_gamma = gamma[index,i]-gamma[index,i+1]
            photon_px += photon_gamma*px[index,i]/np.sqrt(px[index,i]**2+py[index,i]**2)
            photon_py += photon_gamma*py[index,i]/np.sqrt(px[index,i]**2+py[index,i]**2)
        #print(ith,'photon_g=',photon_gamma,'photon_px=',photon_px,'photon_py=',photon_py,'py=',py[index,-1])
        #print(ith,'photon_px=',photon_px,'photon_py=',photon_py,'py=',py[index,-1])
        photon_x[ith]=photon_px;photon_y[ith]=photon_py;
    
    plt.subplot(2,1,1)
    index=np.where(py[:,-1] > 0.0)
    #print(index)
    plt.scatter(photon_y[index],py[index,-1],s=20,c=(192.0/255.0,0.0,0.0),label='py>0',edgecolors='None')
    index=np.where(py[:,-1] <= 0.0)
    plt.scatter(photon_y[index],py[index,-1],s=20,c=(0.0,192.0/255.0,0.0),label='py<=0',edgecolors='None')
    plt.legend(loc='upper right')
    plt.grid()
    #plt.xlim(-2.5,2.5)
    #plt.ylim(-50,50)
    plt.xlabel('photon_$p_y [m_ec]$',fontdict=font)
    plt.ylabel('$p_{yf} [m_ec]$',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    plt.title(r'$\xi_0$='+insert,fontdict=font)
    #plt.show()
    #plt.figure(figsize=(100,100))
 
    plt.subplot(2,1,2)
    plt.scatter(x[:,-1]/2.0/np.pi,y[:,-1]/2.0/np.pi,s=20,c=(192.0/255.0,0.0,0.0),label='QED RR',edgecolors='None')
    plt.legend(loc='upper right')
    #plt.grid()
    plt.xlim(0,100)
    plt.ylim(-50,50)
    plt.xlabel('$x [\mu m]$',fontdict=font)
    plt.ylabel('$y [\mu m]$',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.title(r'$\xi_0$='+insert,fontdict=font)

    fig = plt.gcf()
    fig.set_size_inches(12, 18)
    fig.savefig('./'+insert+'py_vs_ppy.png',format='png',dpi=60)
    plt.close("all")
    #plt.title('radiated total number',fontdict=font)
    #    print('px,',px[index,i-2:i+3],px[index,i+1])
    #    print('py,',py[index,i-2:i+3],py[index,i+1])
    #    print('gamma,',gamma[index,i-2:i+3],gamma[index,i+1])
    #    print(i)
    #(1000.04-945.86)
    #photon_px=



