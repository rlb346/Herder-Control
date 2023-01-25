import numpy as np
from numba import njit
from scipy.interpolate import RectBivariateSpline
import random

import Parameters as P
import Guidance_vector as G

nxy = P.nxy
boundary_scale = P.boundary_scale
length = P.length
t_end = P.t_end
diffusive_timescale = P.diffusive_timescale
dt = P.dt
D = P.D
deltaxy = P.deltaxy
n_stakeholders = P.n_stakeholders
n_particles = P.n_particles
brownian_motion = P.brownian_motion
Dp = P.Dp
mu = P.mu
radius = P.particle_radius


#some preliminary stuff
grid = np.zeros((nxy,nxy))     # solution array
#Because of the way we index variables on the meshgrid, we have to transpose grid when we want to plot it.
#This is because X[0] doesn't give the first x value, but X[:,0] does. Instead of reversing the indices on everything, it is easier to just transpose the final resuult.
Sources = np.zeros((nxy,nxy))  # source term array. 


interval = (boundary_scale-1)/2
x = np.linspace(-interval*length,(1+interval)*length,nxy)
y = np.linspace(-interval*length,(1+interval)*length,nxy)
X,Y = np.meshgrid(x,y) 

nTauRun = t_end/diffusive_timescale              # number of intrinsic timescale to run for
t_end = nTauRun*diffusive_timescale              # Change the end time so it is a multiple of the diffusive timescale.
#will this cause problems with the controller times not lining up? Might have to change this.
nt   = int(np.ceil(t_end/dt)) +1 # number of timesteps
tlarge = np.linspace(0,t_end,nt)
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
def process(uchem,stake_target,grid,time, initialpos,chase_index): 
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
            Sources[xindex[n],yindex[n]] = uchem[n] 
        finite_difference(grid,Sources)
        interp = RectBivariateSpline(x,y,grid) #default for RectBivariateSpline is cubic interpolation
        for m in range(n_particles):
            if m < n_stakeholders: #chase_index might mess it up when chase_index == 0
                #v_desired = guidance_vector(xy[k,:,0],stake_target, xy[k,:,chase_index]) #this is for chased particle, then will need to do for others. 
                v_desired = G.guidance_vector(xy[k,:,0],stake_target, xy[k,:,n_stakeholders:],chase_index)
                vxy[k,0,m] = v_desired[0] + brownian_motion*np.random.normal(0,1)*np.sqrt(2*Dp)/np.sqrt(dt)
                vxy[k,1,m] = v_desired[1] + brownian_motion*np.random.normal(0,1)*np.sqrt(2*Dp)/np.sqrt(dt)     
            else:
                vxy[k,0,m] = mu*(interp(xy[k-1,0,m]+deltaxy,xy[k-1,1,m]) - interp(xy[k-1,0,m]-deltaxy,xy[k-1,1,m]))/(2*deltaxy) + brownian_motion*np.random.normal(0,1)*np.sqrt(2*Dp)/np.sqrt(dt)
                vxy[k,1,m] = mu*(interp(xy[k-1,0,m],xy[k-1,1,m]+deltaxy) - interp(xy[k-1,0,m],xy[k-1,1,m]-deltaxy))/(2*deltaxy) + brownian_motion*np.random.normal(0,1)*np.sqrt(2*Dp)/np.sqrt(dt)
            xy[k+1,:,m] = xy[k,:,m] + vxy[k,:,m]*dt
            #now it just moves particle m, and that seems to work. Fix the hard spheres to do the same.
            #now I'm getting the solver not working sometimes. Why?
            
            #hard sphere interactions.
            xcheck = xy[k+1,0,:] #notice this is not making a copy, it is just attaching a name to this section of the data.
            ycheck = xy[k+1,1,:]
            random.shuffle(particle_list)
            for ii in particle_list:
                if ii != m:
                    dij = np.sqrt((xcheck[ii]-xcheck[m])**2+(ycheck[ii]-ycheck[m])**2)           
                    if dij < 2*radius:
                        xcheck[m] = xcheck[m] - 1*(dij-2*radius)*(xcheck[m]-xcheck[ii])/dij
                        ycheck[m] = ycheck[m] - 1*(dij-2*radius)*(ycheck[m]-ycheck[ii])/dij
    return (xy[-1])
