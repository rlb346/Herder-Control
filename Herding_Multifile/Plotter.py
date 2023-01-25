import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from IPython.display import display, clear_output
import matplotlib.patches as patches

import Parameters as P
import Process as Pro

length = P.length
radius = P.particle_radius
target_type = P.target_type
n_particles = P.n_particles
n_stakeholders = P.n_stakeholders
target_position = P.target_position
folder = P.folder
filename = P.filename
X = Pro.X
Y = Pro.Y
grid = Pro.grid

cm = 1E6 #this is to show the plot in micrometers (used to be centimeters)
Clow = 0
Chigh = 0.1 #0.6
Elow = 0
Ehigh = 0.03
margin = length/10
lower = 0-margin
higher = length+margin
ntail = 30
marksize = 0.75*0.775*2*radius*234/(higher-lower)#72 points per inch times 3.25 inches across=234. 0.775 is to account for axis label 
#the 0.75 is a fudge factor to make it look right. I can't figure out why, but without it the plot is wrong. 

#plotting function
plt.rc("font",size = 10)
colors = ['tab:orange','tab:red','tab:green','tab:olive','tab:brown','tab:cyan','tab:purple'] 


# targett = np.zeros((len(t),2))
# for i in range(len(t)):
#     targett[i] = target_trajectory(t[i])[n_stakeholders:,1] #this might be wrong

def plotter(k,xanswer,yanswer,uresult):
    plt.figure(figsize=(3.25, 3.25), dpi= 200, facecolor='w', edgecolor='k')
    plt.xlabel("Position ($\mu$m)", fontsize = 12) 
    #plt.plot(targett[:380,0]*cm,targett[:380,1]*cm, "r--", lw = 0.75) #why 380? this should be a parameter.
    plt.plot(uresult[0]*cm,uresult[1]*cm,'s',color = 'm', ms = 2)
    if target_type == 'position':
        for c in range(n_particles):
# =============================================================================
#             if c < n_stakeholders:
#                 order = 3
#             else:
#                 order = 2
# =============================================================================
            plt.plot(target_position[c,0]*cm,target_position[c,1]*cm,'x',color = colors[c]) #double check this
            plt.plot(xanswer[k,c]*cm,yanswer[k,c]*cm, "o", color = colors[c], markersize = marksize)
            if k< ntail:
                plt.plot(xanswer[:k+1,c]*cm,yanswer[:k+1,c]*cm,color = colors[c], linewidth = 1.5)
            else:
                plt.plot(xanswer[k-ntail:k+1,c]*cm,yanswer[k-ntail:k+1,c]*cm,color = colors[c], linewidth = 1.5)
        #plt.plot(target_position[:,0]*cm,target_position[:,1]*cm,'x',color = 'r') #double check this
    
    #plot passive particles    
    # plt.plot(xanswer[k,n_stakeholders:]*cm,yanswer[k,n_stakeholders:]*cm, "o", color = "r", markersize = marksize)
    # if k< ntail:
    #     plt.plot(xanswer[:k+1,n_stakeholders:]*cm,yanswer[:k+1,n_stakeholders:]*cm,color = "r", linewidth = 1.5)
    # else:
    #     plt.plot(xanswer[k-ntail:k+1,n_stakeholders:]*cm,yanswer[k-ntail:k+1,n_stakeholders:]*cm,color = "r", linewidth = 1.5)
    # #plt.plot(x0[:,0]*cm, x0[:,1]*cm, "s", markersize = 10,color = "royalblue")
    # #plot stakeholders
    # plt.plot(xanswer[k,:n_stakeholders]*cm,yanswer[k,:n_stakeholders]*cm, "o", color = "tab:orange", markersize = marksize)
    # if k< ntail:
    #     plt.plot(xanswer[:k+1,:n_stakeholders]*cm,yanswer[:k+1,:n_stakeholders]*cm,color = "tab:orange", linewidth = 1.5)
    # else:
    #     plt.plot(xanswer[k-ntail:k+1,:n_stakeholders]*cm,yanswer[k-ntail:k+1,:n_stakeholders]*cm,color = "tab:orange", linewidth = 1.5)
    axes = plt.gca()
    divider = make_axes_locatable(axes)
    # if mode == "chemical":
    plt.contourf(X*cm,Y*cm,grid.T/1000,levels=np.linspace(Clow,Chigh,26), extend = "max") #recall that grid must be transposed to give the correct plot.
    cax = divider.append_axes("left", size = "5%", pad = 0.05)
    colorbar = plt.colorbar(cax = cax, shrink = 0.81, pad = 0.03)
    colorbar.set_label("Concentration ($\dfrac{\mathrm{mol}}{\mathrm{L}}$)", fontsize = 12)   
    cax.yaxis.tick_left()
    cax.yaxis.set_label_position('left')
    axes.set_aspect('equal', adjustable='box')
    axes.set_xlim([lower*cm,higher*cm])
    axes.set_ylim([lower*cm,higher*cm])
    axes.yaxis.tick_right()
    clear_output(wait=True)
    display(plt.gcf())
    kfilled = str(k).zfill(5)
    plt.savefig(f"{folder}/{filename}Plot{kfilled}.png",bbox_inches="tight") 
    plt.close()
    return