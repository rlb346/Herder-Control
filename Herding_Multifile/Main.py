import numpy as np
from scipy.optimize import least_squares
#import importlib
import Parameters as P #Write the name of your parameters file here. 
import Cost_function as C
import Process as Pro
import Model as M
import Plotter
import Save_data as S
#importlib.reload(P) #Reload to see any changes made to the parameter file.


##Parameters
n_times = P.n_times
n_particles = P.n_particles
initial_position = P.initial_position
n_vision = P.n_vision
n_stakeholders = P.n_stakeholders
nxy = P.nxy
t = P.t
target_position = P.target_position
radius = P.particle_radius
radius_h = P.herder_radius
mu = P.mu
bounds = P.bounds
decision_times = P.decision_times
umax = P.umax
dtratio = P.dtratio
save = P.save

tlarge = Pro.tlarge
ncutoffsteps = M.ncutoffsteps
grid = Pro.grid





##run simulation and controller
xy = np.zeros((n_times,2,n_particles))
xy[0] = initial_position.T
vxy = np.zeros((n_times,2,n_particles))
##vxy[0] is zero
uguess = np.zeros(2*(n_vision+1))#+length/2 #oh no, the guess value matters!

ustore = np.zeros((n_times,n_stakeholders))
staketargetstore = np.zeros((n_times,2))
modelpos = xy.copy()
timer = np.zeros(len(t))
switch_targets = True
i = 0
while True: #does python have do while loops? Or move the += to end of this block.
    #switching rule
    distance = np.sqrt((xy[i,0,:]-target_position[:,0])**2 + (xy[i,1,:]-target_position[:,1])**2)
    if (distance[n_stakeholders:] > 2*radius).any(): #2 should be a parameter 
        if switch_targets:
            chase_index = np.argmax(distance[n_stakeholders:])+n_stakeholders
            switch_targets = False
        elif distance[chase_index] < radius:
            switch_targets = True     
            #should I have a continue statement here or something?
    else:
        chase_index = 0
    #heuristic for choosing guess value.
    right = xy[i,0,chase_index]-target_position[chase_index,0] #particle is to the right of target if this is positive
    above = xy[i,1,chase_index]-target_position[chase_index,1] #particle is above target if this is positive
    dist = np.sqrt(right**2 + above**2)
    uguess[::2] = xy[i,0,chase_index]-right/dist*8*radius*np.sign(mu)
    uguess[1::2]= xy[i,1,chase_index]-above/dist*8*radius*np.sign(mu)
    if chase_index == 0 and distance[0] < radius_h: 
        break
    
    stakehistory =  xy[:,:,:n_stakeholders] #the past positions of the stakeholders. Goes into model. Needs a cutoff. 
    if i == 0: #this is to make array sizes work 
        result = least_squares(C.objective,uguess, bounds = bounds,args = (decision_times[:i+n_vision+1],np.zeros((n_stakeholders,0)), xy[i],i,stakehistory[:i],chase_index),verbose = 1) 
    elif i < ncutoffsteps: 
        result = least_squares(C.objective,uguess, bounds = bounds,args = (decision_times[:i+n_vision+1],ustore[:i].T, xy[i],i,stakehistory[:i],chase_index),verbose = 1)
    else:
        result = least_squares(C.objective,uguess, bounds = bounds, args = (decision_times[i-ncutoffsteps:i+n_vision+1],ustore[i-ncutoffsteps:i].T, xy[i],i,stakehistory[i-ncutoffsteps:i],chase_index),verbose = 1)

    uresult = result.x #this gives the input divided by uemax. 
    uresultchem = np.ones((n_stakeholders,n_vision))*umax #this is now in two places. Fix it.
    xy[i+1] = Pro.process(uresultchem[:,0],uresult[:2],grid,tlarge[i*dtratio:(i+1)*dtratio],xy[i],chase_index)
    ustore[i] =  uresultchem[:,0] #this is set up for later, but right now it's awkward.
    staketargetstore[i+1] = uresult[:2]
    
    Plotter.plotter(i+1,xy[:,0],xy[:,1],uresult[:2])
    if save:
        S.save_data(i) #this function is not yet implemented.
    #print(i, result.nfev)
    i += 1


#%%
    print(np.linalg.norm(xy[i,:,0]-staketargetstore[i])/radius)