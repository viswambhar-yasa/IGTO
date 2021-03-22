#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#------------------------------------------------------------------------------------------------------------------#
#PLOTTING - python code used to generate graphs for values obtained after topology optimization and  time analysis.
#------------------------------------------------------------------------------------------------------------------#

import matplotlib.pyplot as plt

import numpy as np
import os


def Folder(path):
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except OSError:
        print('Error: Creating directory. ' + path)
      
def plotting(ii,CC,VV,Mnd,element_density,optimizer,option):
    '''
    This function generates plots from the data obtained after toplogy optimization.

    INPUTS :
    -------
        ii      :Arrat
                 Contains the iteration numbers

        CC       :Array
                 Compliance at each iteration

        VV       :Array
                 Volume fraction at each iteration

        Mnd     : int 
                    Measure of discretness

        element_density : array
                            element density of each element obtained after optimization

        optimizer  : str
                        Name of the optimizer used in toplogy optimization 

        option :int 
                    Boundary condition used Ex-0 (Cantilever beam with load along the bottom edge of the free end) used in saving plots with that name


    Returns 
    ---------
            3 graphs iteration VS compliance , iteratiion VS volume fraction and Measure of discretness(histogram of element density)


    '''

    fig, ax = plt.subplots()
    ax.plot(ii,CC)
    ax.plot(ii[-1], CC[-1],'bo')
    ax.text(ii[-1], CC[-1],str(round(CC[-1],2)),color='b')
    ax.set_xlabel('iteration')
    ax.set_ylabel('Compliance')
    ax.set_xlim(0,ii[-1]+2)
    ax.set_ylim(0,CC[3])
    ax.set_title('Compliance change in each iteration')
    Folder('./results/')
    path="./results/"+optimizer+"_ComplianceVSiteration.png"
    fig.savefig(path)

    fig, ax1 = plt.subplots()
    ax1.plot(ii,VV)
    ax1.plot(ii[-1], VV[-1],'bo')
    ax1.text(ii[-1], VV[-1],str(round(VV[-1],2)),color='b')
    ax1.set_xlabel('iteration') 
    ax1.set_ylabel('Volume fraction')
    ax1.set_xlim(0,ii[-1]+2)
#ax.set_ylim(volume_frac-0.1,volume_frac+0.1)
    ax1.set_title('Volume fraction change in each iteration')
    Folder('./results/')
    path="./results/"+optimizer+"_Volume_fractionVSiteration.png"
    fig.savefig(path)
    
    fig, ax2 = plt.subplots()
    number,bins,pat=ax2.hist(element_density,bins=10)
    
    ax2.text(0.25,max(number)+0.5,'Measure of discreteness: '+str(round(Mnd,2)))
    ax2.set_xlabel('volume fraction ') 
    ax2.set_ylabel('Number of elements')
#ax.set_ylim(volume_frac-0.1,volume_frac+0.1)
    ax2.set_title('Discretness of volume fraction')
    Folder('./results/')
    path="./results/"+optimizer+"_discretness.png"
    fig.savefig(path)


def time_analysis_plotting(x,y,optimizer):
    fig, ax = plt.subplots() 
    labels = x
    data = y
    nn=len(data)

    explode=np.ones(nn)*0.02
    
    plt.pie(data, labels=labels, autopct='%1.1f%%', startangle=90, pctdistance=0.85, explode = explode)
    #draw circle
    circle = plt.Circle((0,0),0.70,fc='white')
    fig = plt.gcf()
    fig.gca().add_artist(circle)
    # Equal aspect ratio ensures that pie is drawn as a circle
    ax.axis('equal')  
    path="./results/"+optimizer+"_time_analysis.png"
    ax.set_title("Time analysis Function calls") 

    fig.savefig(path)
    pass






















