import numpy as np
import Parameters as P #Write the name of your parameters file here. 
import Process as Pro
import Plotter
import Save_data as S
import Guidance_vector as G
import time

import pandas as pd


##Parameters
n_particles = P.n_particles
initial_position = P.initial_position
n_herders = P.n_herders
nxy = P.nxy
t = P.t
target_position = P.target_position
radius = P.particle_radius
radius_h = P.herder_radius
mu = P.mu
umax = P.umax
dtratio = P.dtratio
save = P.save
d_tol = P.d_tol
d_prec = P.d_prec
R_close = P.R_close
delta_t = P.delta_t

tsmall = Pro.tsmall
grid = Pro.grid



## Run simulation and controller
xy = np.array([initial_position.T])
switch_targets = True
finished = False
chase_index = 1
i = 0
storeerrors = []
finished_counter = 0
while True: # Run the loop until it triggers the break statement.
    target_position = P.matchTargets(target_position)
    # switching rule
    time0 = time.perf_counter()
    distance = np.sqrt((xy[i,0,:] - target_position[:,0])**2 + (xy[i,1,:] - target_position[:,1])**2)
    if (distance[n_herders:] > d_tol).any(): #       
        nothing = 0 #This is a useless line of code
    else:
        if finished == False:
            print("Completed")
            finished = True   
    if finished:
        finished_counter += 1
    if finished_counter == 3000: #run for 5 more minutes after finishing
        print("5 minutes completed")
        break
        
    if distance[chase_index] < d_prec:
        switch_targets = True 
    if switch_targets:
        chase_index = np.argmax(distance[n_herders:])+n_herders
        switch_targets = False   
# =============================================================================
#     if distance[0] < 1.5*(radius+radius_h) and finished: 
#         print("Completed")
#         break
# =============================================================================
    if i >= 40000:
        print("Too many iterations")
        break
    
    
    #heuristic for choosing herder position.
    error = xy[i,:,chase_index] - target_position[chase_index,:]
    dist = np.linalg.norm(error)
    herder_target = xy[i,:,chase_index] - error/dist*(R_close + radius + radius_h)*np.sign(mu)
    
    time1 = time.perf_counter()
    herder_velocity = G.guidance_vector(xy[i,:,:n_herders].flatten(), herder_target, xy[i,:,n_herders:], chase_index)
    
    time2 = time.perf_counter()
    xnew = [Pro.process(umax, herder_velocity, grid, t[i] + tsmall, xy[i])]
    time3 = time.perf_counter()
    xy = np.concatenate((xy, xnew))
    t = np.append(t,t[-1] + delta_t)
    
    
    if i%20 == 0:
        Plotter.plotter(i, xy[:,0], xy[:,1], herder_target)
        errorsall = xy[i,:,n_herders:] - target_position[n_herders:,:].T
        sumerrors = np.linalg.norm(errorsall/radius)
        print(i,np.linalg.norm(xy[i,:,0] - herder_target)/radius,sumerrors)
        #print(time1-time0,time2-time1,time3-time2)
        storeerrors.append(sumerrors)
        if save:
            S.save_data(i,t,xy) 
    i += 1
    


#%%
print(np.amax(grid)/1000)
