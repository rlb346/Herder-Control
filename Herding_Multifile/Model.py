import numpy as np
from numba import njit
from scipy.optimize import fsolve

import Parameters as P

cutoff_percent = P.cutoff_percent
length = P.length
D = P.D
delta_t = P.delta_t
n_stakeholders = P.n_stakeholders
radius = P.particle_radius
n_particles = P.n_particles
t = P.t


@njit #The chemical kernel for a single step function, or in other words, one half of delta h. 
def f(t,tn,r,u_i,mu, D): 
    if t>tn:
        rsquared = r**2
        return -mu*u_i*np.exp(-rsquared/(4*D*(t-tn)))/(2*D*np.pi*rsquared) #note that this is missing the (x-x0), added in the chemical_velocity function.
    else:
        return 0

# Find the cutoff time.   
def find_cutoff(tdiff, delta_t):
    def h0_derivative(t): #The derivative of delta h.
        if t <= delta_t:
            return 1E8 #Set to a big number so the solver avoids it, this is how I implement the Heaviside function.
        else:
            return np.exp(-tdiff/t)/t**2-np.exp(-tdiff/(t-delta_t))/(t-delta_t)**2
    maximum = fsolve(h0_derivative,0.5*tdiff + 3/4*delta_t) 
    def h0(t):
        return np.exp(-tdiff/t)-np.exp(-tdiff/(t-delta_t))
    def function_to_solve(t):
        return h0(t) - cutoff_percent*h0(maximum)
    answer = fsolve(function_to_solve,8*tdiff) #If it doesn't work, try a new guess value here
    return answer[0]
tdiff = length**2/4/D
tcutoff = find_cutoff(4*tdiff,delta_t)
if tcutoff < 0:
    raise     #If this raises an exception, try a new guess value above
ncutoffsteps = int(tcutoff/delta_t) #number of timesteps before cutoff. Note that this rounds down. Maybe you don't want that.



@njit #finds chemical contribution to vx and vy
def chemical_velocity(xy,t,stakepos,u,decision_times, mu, D,stakehistory): 
    x = xy[0]
    y = xy[1]
    vx = 0.0
    vy = 0.0
    for i in range(len(stakehistory)+1): 
        for j in range(n_stakeholders):
            if i == len(stakehistory):#Don't actually need this because decision_times>t here. #It actually does get used, what did I mean by that?
                rad = np.sqrt((stakepos[0,j]-x)**2 + (stakepos[1,j]-y)**2)
                if rad < radius:
                    continue
                vx += f(t,decision_times[i],rad,u[j,i],mu,D)*(x-stakepos[0,j])
                vx -= f(t,decision_times[i] + delta_t,rad,u[j,i],mu,D)*(x-stakepos[0,j])
                vy += f(t,decision_times[i],rad,u[j,i],mu,D)*(y-stakepos[1,j]) 
                vy -= f(t,decision_times[i] + delta_t,rad,u[j,i],mu,D)*(y-stakepos[1,j])
            else:  
                rad = np.sqrt((stakehistory[i,0,j]-x)**2 +(stakehistory[i,1,j]-y)**2)
                if rad < radius:
                    continue
                vx += f(t,decision_times[i],rad,u[j,i],mu,D)*(x-stakehistory[i,0,j])
                vx -= f(t,decision_times[i] + delta_t,rad,u[j,i],mu,D)*(x-stakehistory[i,0,j])
                vy += f(t,decision_times[i],rad,u[j,i],mu,D)*(y-stakehistory[i,1,j])
                vy -= f(t,decision_times[i] + delta_t,rad,u[j,i],mu,D)*(y-stakehistory[i,1,j])
    return np.array([vx,vy])

@njit #model for particle motion given a set of reaction rates and voltages
def model(uchem,xstake,decision_times, measurement,start, stop, mu, D, stakehistory): 
    xy  = np.zeros((stop-start + 1,2,n_particles))
    dxy = np.zeros((stop-start + 1,2,n_particles)) 
    xy[0] = measurement
    #initial time point
    for j in range(n_particles):
        dxy[0,:,j]  = chemical_velocity(xy[0,:,j],t[start],xy[0,:,:n_stakeholders],uchem,decision_times,mu, D,stakehistory) #maybe speed it up by only passing a slice of u
    #loop over time 
    stakehistoryexpanded = stakehistory
    for i in range(1,stop-start+1): 
        xy[i,:,n_stakeholders:] = xy[i-1,:,n_stakeholders:] + dxy[i-1,:,n_stakeholders:]*delta_t #Explicit Euler                 
        stakehistoryexpanded = np.append(stakehistoryexpanded,np.expand_dims(xy[i-1,:,:n_stakeholders],0),0)
        for j in range(n_particles):
            if j < n_stakeholders:
                #dxy[i,:,j] = dstake[i]
                xy[i,:,j] = xstake[i-1] #not built for j > 1
            else:
                dxy[i,:,j]  = chemical_velocity(xy[i,:,j],t[i+start],xy[i,:,:n_stakeholders],uchem,decision_times,mu, D,stakehistoryexpanded)    
    return (dxy,xy)
