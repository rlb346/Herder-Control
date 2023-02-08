import numpy as np
from numba import njit
from scipy.interpolate import RectBivariateSpline
import random

import Parameters as P


nxy = P.nxy
boundary_scale = P.boundary_scale
length = P.length
#t_end = P.t_end
#diffusive_timescale = P.diffusive_timescale
dt = P.dt
D = P.D
deltaxy = P.deltaxy
n_stakeholders = P.n_stakeholders
n_particles = P.n_particles
brownian_motion = P.brownian_motion
Dp = P.Dp
Dp_h = P.Dp_herder
mu = P.mu
radius = P.particle_radius
radius_h = P.herder_radius
dtratio = P.dtratio
delta_t = P.delta_t


#some preliminary stuff
grid = np.zeros((nxy,nxy))     # solution array
#Because of the way we index variables on the meshgrid, we have to transpose grid when we want to plot it.
#This is because X[0] doesn't give the first x value, but X[:,0] does. Instead of reversing the indices on everything, it is easier to just transpose the final resuult.
Sources = np.zeros((nxy,nxy))  # source term array. 


interval = (boundary_scale-1)/2
x = np.linspace(-interval*length,(1+interval)*length,nxy)
y = np.linspace(-interval*length,(1+interval)*length,nxy)
X,Y = np.meshgrid(x,y) 

#nTauRun = t_end/diffusive_timescale              # number of intrinsic timescale to run for
#t_end = nTauRun*diffusive_timescale              # Change the end time so it is a multiple of the diffusive timescale.
#will this cause problems with the controller times not lining up? Might have to change this.
#nt   = int(np.ceil(t_end/dt)) +1 # number of timesteps
#tlarge = np.linspace(0,t_end,nt)

tsmall = np.linspace(0,delta_t,dtratio)
particle_list = list(range(n_particles))



@njit #this does one iteration of the Gauss-Seidel algorithm 
def finite_difference(grid,Sources):
    for iindex in range(1,nxy-1):
        for jindex in range(1,nxy-1):
            grid[iindex,jindex] += (
                D*dt/deltaxy**2*(grid[iindex-1,jindex] - 2*grid[iindex,jindex] + grid[iindex+1,jindex]) + 
                D*dt/deltaxy**2*(grid[iindex,jindex-1] - 2*grid[iindex,jindex] + grid[iindex,jindex+1]) + 
                dt*Sources[iindex,jindex]/deltaxy**2 )#Divide the last term by deltaxy^2 to get concentration.
              #Zero boundary conditions
    return

xindex = np.zeros(n_stakeholders, dtype = int)
yindex = np.zeros(n_stakeholders, dtype = int)

#this function runs the physics and gives the resulting positions
def process(uchem,stake_velocity,grid,time, initialpos): 
    xy = np.zeros((len(time),2,n_particles))
    xy[0] = initialpos
    vxy = np.zeros((len(time),2,n_particles))
    for k in range(len(time)-1): 
        for n in range(n_stakeholders):
            #Set previous stakeholder source term back to zero.
            Sources[xindex[n],yindex[n]] = 0  
            #Find the new location for the stakeholder source term.
            xindex[n] = np.argmin((X[0,:]-xy[k,0,n])**2)
            yindex[n] = np.argmin((Y[:,0]-xy[k,1,n])**2)
            Sources[xindex[n],yindex[n]] = uchem#[n] 
        finite_difference(grid,Sources)
        interp = RectBivariateSpline(x,y,grid) #default for RectBivariateSpline is cubic interpolation
        for m in range(n_particles): #does this cause errors by moving them in order?
            if m < n_stakeholders: 
                v_desired = stake_velocity
                vxy[k,0,m] = v_desired[0] + brownian_motion*np.random.normal(0,1)*np.sqrt(2*Dp)/np.sqrt(dt)
                vxy[k,1,m] = v_desired[1] + brownian_motion*np.random.normal(0,1)*np.sqrt(2*Dp)/np.sqrt(dt)     
            else:
                vxy[k,0,m] = mu*(interp(xy[k-1,0,m]+deltaxy,xy[k-1,1,m]) - interp(xy[k-1,0,m]-deltaxy,xy[k-1,1,m]))/(2*deltaxy) + brownian_motion*np.random.normal(0,1)*np.sqrt(2*Dp_h)/np.sqrt(dt)
                vxy[k,1,m] = mu*(interp(xy[k-1,0,m],xy[k-1,1,m]+deltaxy) - interp(xy[k-1,0,m],xy[k-1,1,m]-deltaxy))/(2*deltaxy) + brownian_motion*np.random.normal(0,1)*np.sqrt(2*Dp_h)/np.sqrt(dt)
            xy[k+1,:,m] = xy[k,:,m] + vxy[k,:,m]*dt
            
            #hard sphere interactions.
            xcheck = xy[k+1,0,:] #notice this is not making a copy, it is just attaching a name to this section of the data.
            ycheck = xy[k+1,1,:]
            #random.shuffle(particle_list) #hard to shuffle with different radius.
            #implement this as seperate numba function.
            for ii in particle_list:
                if ii != m:
                    dij = np.sqrt((xcheck[ii]-xcheck[m])**2+(ycheck[ii]-ycheck[m])**2) 
                    if ii < n_stakeholders:
                        overlap_distance = radius+radius_h
                    else:
                        overlap_distance = 2*radius
                    if dij < overlap_distance:
                        xcheck[m] = xcheck[m] - 1*(dij-overlap_distance)*(xcheck[m]-xcheck[ii])/dij
                        ycheck[m] = ycheck[m] - 1*(dij-overlap_distance)*(ycheck[m]-ycheck[ii])/dij
    return (xy[-1])

