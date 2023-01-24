import numpy as np
from numba import njit

import Parameters as P


target_type = P.target_type
length = P.length
n_particles = P.n_particles
feedback_gain = P.feedback_gain

##target
if target_type == 'circle':
    v0 = P.v0
    rd = P.rd
    xt = length/2         #center of the circle
    yt = length/2
    @njit #desired trajectory
    def target_trajectory(t):
        x = np.zeros(n_particles)
        y = np.zeros(n_particles)
        # x[0] = xt+rd*np.cos(v0*t/rd)
        # y[0] = yt+rd*np.sin(v0*t/rd)
        x[1] = xt+rd*np.cos(v0*t/rd) #fix this so it's not zero
        y[1] = yt+rd*np.sin(v0*t/rd)
        return np.vstack((x,y))

    @njit #velocity target (derivative of desired trajectory with a feedback correction)
    def target_velocity(xs,t,xyought):
        k = feedback_gain
        x = xs[0]
        y = xs[1]
        xought = xyought[0]
        yought = xyought[1]
        dx = np.zeros(n_particles)
        dy = np.zeros(n_particles)
        # dx[0] = -v0*np.sin(v0*t/rd)-k*(x[0]-xought[0])
        # dy[0] = v0*np.cos(v0*t/rd) -k*(y[0]-yought[0])
        dx[1] = -v0*np.sin(v0*t/rd) -k*(x[1]-xought[1])
        dy[1] =v0*np.cos(v0*t/rd) -k*(y[1]-yought[1])
        return np.vstack((dx,dy))
elif target_type == 'figure8':
    v0 = P.v0
    rd = P.rd
    xt = length/2         #center of the figure-eight
    yt = length/2
    @njit #desired trajectory
    def target_trajectory(t):
        x = np.zeros(n_particles)
        y = np.zeros(n_particles)
        # x[0] = xt+rd*np.cos(v0*t/rd)
        # y[0] = yt+rd*np.sin(v0*t/rd)
        x[1] = xt+rd*np.sin(v0*t/rd)
        y[1] = yt+rd*np.sin(v0*t/rd)*np.cos(v0*t/rd)
        return np.vstack((x,y))

    @njit #velocity target (derivative of desired trajectory with a feedback correction)
    def target_velocity(xs,t,xyought):
        k = feedback_gain
        x = xs[0]
        y = xs[1]
        xought = xyought[0]
        yought = xyought[1]
        dx = np.zeros(n_particles)
        dy = np.zeros(n_particles)
        # dx[0] = -v0*np.sin(v0*t/rd)-k*(x[0]-xought[0])
        # dy[0] = v0*np.cos(v0*t/rd) -k*(y[0]-yought[0])
        dx[1] = v0*np.cos(v0*t/rd) -k*(x[1]-xought[1])
        dy[1] = v0*(-np.sin(v0*t/rd)**2+np.cos(v0*t/rd)**2) -k*(y[1]-yought[1])
        return np.vstack((dx,dy))
elif target_type == 'position':
    target_position = P.target_position
    # for i in range(1,len(target_position)):
    #     target_position[i] = length/4+np.random.rand(2)*length/2 #for testing
    @njit
    def target_trajectory(t): #check if this is right dimensions
        return target_position
    
# =============================================================================
#     @njit
#     def target_velocity(xs,t,chase_index): 
#         #k = feedback_gain
#         x = xs[0]
#         y = xs[1]    
#         dx = np.zeros(n_particles)
#         dy = np.zeros(n_particles)
#         distance = np.sqrt((x[chase_index]-target_position[chase_index,0])**2 + (y[chase_index]-target_position[chase_index,1])**2)
#         dx[chase_index] = v0*(target_position[chase_index,0]-x[chase_index])/distance
#         dy[chase_index] = v0*(target_position[chase_index,1]-y[chase_index])/distance
#         #the rest are zero
#         return np.vstack((dx,dy))
# =============================================================================
else:
    raise Exception("target_type must be 'circle' or 'figure8'.")