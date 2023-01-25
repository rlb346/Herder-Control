import numpy as np
from numba import njit

import Parameters as P
import Model as M
import Target as T
n_vision = P.n_vision
u_stakeholders = P.n_stakeholders
umax = P.umax
n_stakeholders = P.n_stakeholders
mu = P.mu
D = P.D
n_particles = P.n_particles
t = P.t
k_close = P.k_close
radius = P.particle_radius
W_1 = P.W_1
length = P.length
W_3 = P.W_3

#FIXME get this njit to work 
#@njit #cost function 
def objective(u,decision_times, previousu,xypos,start, stakehistory,chase_index): 
    xstake = np.reshape(u,(n_vision+1,2))
    uchem = np.ones((n_stakeholders,n_vision+1))*umax #in the future this might be calculated, but for now it's static.
    temp = np.hstack((previousu[0,:],uchem[0,:]))
    uchemcombined = np.expand_dims(temp,0)
    for pro in range(1,n_stakeholders):
        temp = np.hstack((previousu[pro,:],uchem[pro,:]))
        temp1 = np.expand_dims(temp,0)
        uchemcombined = np.vstack((uchemcombined,temp1)) 
    vxy,xys = M.model(uchemcombined,xstake,decision_times,xypos,start, start+n_vision, mu, D, stakehistory)  
    
    tooclose = np.zeros(n_vision+1) #tooclose[0] is 0, could change that.
    for k in range(1,n_vision+1): #n_vision into future, does not include initial point.
        rij = np.sqrt((xys[k,0,0] - xys[k,0,chase_index])**2 + (xys[k,1,0] - xys[k,1,chase_index])**2)
        tooclose[k] = max(0,(k_close*2*radius-rij)) 
    residual = xys[:,:,chase_index]-T.target_position[chase_index,:]  #do first and second actually matter?
    return np.concatenate((W_1*residual.flatten()/length,W_3*tooclose.flatten()/radius)) 
