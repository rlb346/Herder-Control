import numpy as np
from scipy.optimize import least_squares
import Parameters as P #Write the name of your parameters file here. 
#import Cost_function as C
import Process as Pro
#import Model as M
import Plotter
import Save_data as S
import Guidance_vector as G


##Parameters
n_particles = P.n_particles
initial_position = P.initial_position
#n_vision = P.n_vision
n_stakeholders = P.n_stakeholders
nxy = P.nxy
t = P.t
target_position = P.target_position
radius = P.particle_radius
radius_h = P.herder_radius
mu = P.mu
#bounds = P.bounds
#decision_times = P.decision_times
umax = P.umax
dtratio = P.dtratio
save = P.save
d_precision = P.d_precision
R_close = P.R_close
delta_t = P.delta_t

tsmall = Pro.tsmall
#ncutoffsteps = M.ncutoffsteps
grid = Pro.grid





##run simulation and controller
xy = np.array([initial_position.T])

#uresultchem = np.ones((n_stakeholders,n_vision))*umax #this is now in two places. Fix it.
#modelpos = xy.copy()
#timer = np.zeros(len(t))
switch_targets = True
chase_index = 0
i = 0
while True: #does python have do while loops? Or move the += to end of this block.
    #switching rule
    distance = np.sqrt((xy[i,0,:]-target_position[:,0])**2 + (xy[i,1,:]-target_position[:,1])**2)
    if (distance[n_stakeholders:] > d_precision).any(): #
        if distance[chase_index] < radius:
            switch_targets = True 
        if switch_targets:
            chase_index = np.argmax(distance[n_stakeholders:])+n_stakeholders
            switch_targets = False   
    else:
        chase_index = 0  
    if distance[0] < 2*radius_h:
        break
    if i >= 30000:
        break
    
    
    #heuristic for choosing herder position.
    error = xy[i,:,chase_index]-target_position[chase_index,:]
    dist = np.linalg.norm(error)
    stake_target = xy[i,:,chase_index] - error/dist*(R_close + radius+radius_h)*np.sign(mu)
    
    stake_velocity = G.guidance_vector(xy[i,:,:n_stakeholders].flatten(),stake_target, xy[i,:,n_stakeholders:],chase_index)
    xnew = [Pro.process(umax,stake_velocity,grid,t[i] + tsmall ,xy[i])]
    #uresultchem[:,0]
    xy = np.concatenate((xy,xnew))
    t = np.append(t,t[-1]+delta_t)
    
    if i%10 == 0:
        print(i,np.linalg.norm(xy[i,:,0]-stake_target)/radius)
        Plotter.plotter(i,xy[:,0],xy[:,1],stake_target)
    if save:
        S.save_data(i) #this function is not yet implemented.
    i += 1
    


#%%
