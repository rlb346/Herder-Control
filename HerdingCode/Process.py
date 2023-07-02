import numpy as np
from numba import njit
from scipy.interpolate import RectBivariateSpline

# Establish parameters 
import Parameters as P
nxy = P.nxy
boundary_scale = P.boundary_scale
length = P.length
dt = P.dt
D = P.D
deltaxy = P.deltaxy
n_herders = P.n_herders
n_particles = P.n_particles
brownian_motion = P.brownian_motion
Dp = P.Dp
Dp_h = P.Dp_herder
mu = P.mu
radius = P.particle_radius
radius_h = P.herder_radius
dtratio = P.dtratio
delta_t = P.delta_t
m_image = P.m_image

#np.random.seed(1)

# Preliminary setup
grid = np.zeros((nxy,nxy))     # solution array
#Because of the way we index variables on the meshgrid, we have to transpose grid when we want to plot it.
#This is because X[0] doesn't give the first x value, but X[:,0] does. Instead of reversing the indices on everything, it is easier to just transpose the final resuult.
Sources = np.zeros((nxy,nxy))  # source term array. 


interval = (boundary_scale-1)/2
x = np.linspace(-interval*length,(1+interval)*length,nxy)
y = np.linspace(-interval*length,(1+interval)*length,nxy)
X,Y = np.meshgrid(x,y) 

tsmall = np.linspace(0,delta_t,dtratio)
particle_list = list(range(n_particles))


#When debugging, comment out the njit.
#@njit #this function does one iteration of the Gauss-Seidel algorithm 
def finite_difference(grid,Sources):
    for iindex in range(1,nxy-1):
        for jindex in range(1,nxy-1):
            grid[iindex,jindex] += (
                D*dt/deltaxy**2*(grid[iindex-1,jindex] - 2*grid[iindex,jindex] + grid[iindex+1,jindex]) + 
                D*dt/deltaxy**2*(grid[iindex,jindex-1] - 2*grid[iindex,jindex] + grid[iindex,jindex+1]) + 
                dt*Sources[iindex,jindex]/deltaxy**2 )#Divide the last term by deltaxy^2 to get concentration.
                #Zero boundary conditions (not written, because they're zero.)
    return

#@njit
def set_sources(Sources, xindex, yindex, xy, k, uchem):
    for n in range(n_herders): 
        #Set previous herder source term back to zero.
        Sources[xindex[n],yindex[n]] = 0  
        #Find the new location for the herder source term.
        xindex[n] = np.argmin((X[0,:] - xy[k,0,n])**2)
        yindex[n] = np.argmin((Y[:,0] - xy[k,1,n])**2)
        Sources[xindex[n],yindex[n]] = uchem 
    return

@njit
def hard_sphere_interactions(particle_list, m, xcheck, ycheck):
    for ii in particle_list:
        if ii != m:
            dij = np.sqrt((xcheck[ii] - xcheck[m])**2 + (ycheck[ii] - ycheck[m])**2) 
            if ii < n_herders:
                overlap_distance = radius + radius_h
            else:
                overlap_distance = 2 * radius
            if dij < overlap_distance:
                xcheck[m] = xcheck[m] - 1*(dij-overlap_distance)*(xcheck[m]-xcheck[ii])/dij
                ycheck[m] = ycheck[m] - 1*(dij-overlap_distance)*(ycheck[m]-ycheck[ii])/dij
    return

@njit
def findgradient(xy, herderposition, uchem): #give it a follower position and a herder position.
    r = np.linalg.norm(herderposition - xy)
    return m_image * uchem / (4 * np.pi * D) * (herderposition - xy) / r**3

xindex = np.zeros(n_herders, dtype = int)
yindex = np.zeros(n_herders, dtype = int)

#this function runs the physics and gives the resulting positions
def process(uchem, herder_velocity, grid, time, initialpos): 
    xy = np.zeros((len(time),2,n_particles))
    xy[0] = initialpos
    vxy = np.zeros((len(time),2,n_particles))
    for k in range(len(time)-1): 
        #set_sources(Sources,xindex,yindex,xy,k,uchem)
        #finite_difference(grid,Sources)
        #interp = RectBivariateSpline(x,y,grid) #default for RectBivariateSpline is cubic interpolation
        for m in range(n_particles): #does this cause errors by moving them in order?
            if m < n_herders: 
                v_desired = herder_velocity
                vxy[k,0,m] = v_desired[0] + brownian_motion * np.random.normal(0,1) * np.sqrt(2*Dp) / np.sqrt(dt)
                vxy[k,1,m] = v_desired[1] + brownian_motion * np.random.normal(0,1) * np.sqrt(2*Dp) / np.sqrt(dt)     
            else:
                mygradient = findgradient(xy[k,:,m],xy[k,:,0],uchem) #this just does first herder
                #fixme: loop over all herders
                vxy[k,0,m] = mu * mygradient[0] + brownian_motion * np.random.normal(0,1) * np.sqrt(2*Dp_h) / np.sqrt(dt)
                vxy[k,1,m] = mu * mygradient[1] + brownian_motion * np.random.normal(0,1) * np.sqrt(2*Dp_h) / np.sqrt(dt)
            xy[k+1,:,m] = xy[k,:,m] + vxy[k,:,m]*dt
            
            #hard sphere interactions.
            xcheck = xy[k+1,0,:] #notice this is not making a copy, it is just attaching a name to this section of the data.
            ycheck = xy[k+1,1,:]
            np.random.shuffle(particle_list) #shuffle so hard sphere interactions are in random order      
            hard_sphere_interactions(particle_list,m,xcheck,ycheck)

    return (xy[-1])

