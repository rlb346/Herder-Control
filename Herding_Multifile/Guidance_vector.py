import numpy as np

import Parameters as P

length = P.length
v0_stakeholder = P.v0_stakeholder
radius = P.particle_radius

little_r = P.little_r 
big_R = P.big_R    
G_path = P.G_path
H_path = P.H_path
G_obstacle = P.G_obstacle
H_obstacle = P.H_obstacle    
#This function gives a guidance vector field to move the stakeholder to the target while avoiding the obstacles.
#The paper it is based on moves the vehicle to the origin, so you must shift the coordinates so the target is the origin.
#@njit
def guidance_vector(position,target,obstacle,chase_index): #the position of the stakeholder and the point it is trying to reach and the location of the obstacle
    x = position[0]/length-target[0]/length #shift coordinate frame so target is the origin, and scale by length
    y = position[1]/length-target[1]/length
    xc = np.zeros(len(obstacle[0]))
    yc = xc.copy()
    xbar = xc.copy()
    ybar = xc.copy()
    d = xc.copy()
    P = xc.copy()
    
    for i in range(len(obstacle[0])): #this might break when you get back to a single passive particle. Oh well. 
        xc[i] = obstacle[0,i]/length - target[0]/length 
        yc[i] = obstacle[1,i]/length - target[1]/length
        xbar[i] = x-xc[i]
        ybar[i] = y-yc[i]
        d[i] = np.sqrt(xbar[i]**2+ybar[i]**2)
        P[i] = 1-np.tanh(2*np.pi*d[i]/(big_R/length)-np.pi)
    delta = np.arctan2(y,x)
    V_circ = np.array([np.sin(delta),-np.cos(delta)])
    V_conv = -1/np.sqrt((np.cos(delta)*x+np.sin(delta)*y)**2)*np.array([x*np.cos(delta)**2+np.cos(delta)*np.sin(delta)*y,y*np.sin(delta)**2+np.cos(delta)*np.sin(delta)*x])
    V_path = G_path*V_conv + H_path*V_circ
    Vg = V_path
    for i in range(len(obstacle[0])):
        #if i == chase_index: #check if this works. If it does, rewrite without loop.
        cdist = np.sqrt(xc[i]**2 + yc[i]**2)
        if cdist > 2*radius/length or i == chase_index:
            Vo_conv = -1/np.sqrt(xbar[i]**4+ybar[i]**4+2*xbar[i]**2*ybar[i]**2-2*little_r**2*xbar[i]**2-2*little_r**2*ybar[i]**2+little_r**2)*np.array([2*xbar[i]**3+2*xbar[i]*ybar[i]**2-2*little_r**2*xbar[i],2*ybar[i]**3+2*xbar[i]**2*ybar[i]-2*little_r**2*ybar[i]])
            Vo_circ = np.array([2*(y-yc[i]),2*(xc[i]-x)])
            V_obstacle = G_obstacle*Vo_conv+H_obstacle*Vo_circ
            Vg += P[i]*V_obstacle
    velocity = Vg
    return velocity/np.linalg.norm(velocity)*v0_stakeholder #try normalizing it, see if it works.    
