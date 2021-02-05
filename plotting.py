import matplotlib.pyplot as plt
from topology import ii,CC,VV

import os


def Folder(path):
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except OSError:
        print('Error: Creating directory. ' + path)
def plotting(ii,CC,element_density,optimizer,option):
    fig, ax = plt.subplots()
    ax.plot(ii,CC)
    ax.set_xlabel('iteration')
    ax.set_ylabel('Compliance')
    ax.set_xlim(0,ii[-1]+2)
#ax.set_ylim(0,CC)
    ax.set_title('Compliance change in each iteration')
    Folder('./results/')
    path="./results/"+optimizer+"_ComplianceVSiteration.png"
    fig.savefig(path)

    fig, ax1 = plt.subplots()
    ax1.plot(ii,VV)
    ax1.set_xlabel('iteration') 
    ax1.set_ylabel('Volume fraction')
    ax1.set_xlim(0,ii[-1]+2)
#ax.set_ylim(volume_frac-0.1,volume_frac+0.1)
    ax1.set_title('Volume fraction change in each iteration')
    Folder('./results/')
    path="./results/"+optimizer+"_Volume_fractionVSiteration.png"
    fig.savefig(path)
    
    fig, ax2 = plt.subplots()
    ax2.hist(element_density,bins=10)
    ax2.set_xlabel('volume fraction ') 
    ax2.set_ylabel('Number of elements')
#ax.set_ylim(volume_frac-0.1,volume_frac+0.1)
    ax2.set_title('Discretness of volume fraction')
    Folder('./results/')
    path="./results/"+optimizer+"_discretness.png"
    fig.savefig(path)



























