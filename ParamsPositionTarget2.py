import numpy as np

##Parameters
folder = 'Data2'
save = False #save the data files or not
filename = 'PositionTarget'

#simulation parameters
n_times = 3001      #number of simulation time steps
t_end = 3000        #end time in seconds
length = 0.0002   #domain length in meters. 
n_electrodes = 0   #number of electrodes
boundary_scale = 3 #number of domain lengths extended for finite difference. 3 is one domain and another on each side
dt = 0.0008       #timestep for physics

#target function parameters
target_type = 'position' #note that the default start position for circle is right side, while default for figure8 is center.
v0 = 1E-6                               #loiter velocity, m/s
v0_stakeholder = 20*v0 #this is supposed to be the speed the stakeholder moves at
xt = length/2                           #center of the circle
yt = length/2
rd = 3.0E-5                             #desired loiter distance, m
# initial_position = np.array([[xt+2.0*rd,yt+1.0*rd],[xt+rd,yt+2*rd],[xt-rd,yt-rd],[xt+rd,yt-2*rd],[xt+rd,yt],[0,0],[0,0],[0,0],[0,0]]) #initial positions of particles
initial_position = np.zeros((12,2))
target_position = np.array([[0,0],[0.0001,0.0001],[0.0001,0.0001],[0.0001,0.0001],[0.0001,0.0001],[0.0001,0.0001],[0.0001,0.0001],[0.0001,0.0001],[0.0001,0.0001]])
target_position = np.array([[0,0]])
for i in range(len(initial_position)-1):
    target_position = np.vstack((target_position,np.array([[0.0001,0.0001]])))
n_stakeholders = 1                       #number of particles designated as stakeholders. Stakeholders are the first elements of initial_position.

#physical parameters
mu = -2.0E-10/1000                    #diffusiophoretic mobility of particle m^2/M s. Divide by 1000 to turn Liters to meters^3.                           
mu_e_stakeholder = 2*2.0E-8          #electrophoretic mobility of particles in Coulomb/Newton * m/s or m^2/V s
epsilon = 78.4*8.8542E-12            #Dielectric constant of water in C^2/Nm^2 
D = 2.30E-9                          #diffusion coefficient of solute. This number is the self-diffusion coeff of water at 25 C.
D = 1.38E-9 #urea in water at 25C. 
brownian_motion = 1                            #Brownian motion. Set to 1 for on or 0 for off
T = 298                              #Temperature in K
viscosity = 8.9E-4                   #viscosity of water
particle_radius = 3E-6                        #m, radius of colloidal particle

#controller parameters
n_vision = 3#was 5         #prediction horizon
W_1 = 20             #weight for residuals in objective function. This is actually W1^2, the way it is written in the paper.
W_2 = 1
W_3 = 1000
cutoff_percent = 0.2 #higher value will solve faster but be less accurate.
feedback_gain = 0.1  #gain for target function
k_close = 1.1
upper_bound_factor = 50       #multiple of uemax to set upper bound at. 1 is normal for a single particle.
strength = 0.5*2*particle_radius/length #I think this would be the strength needed for the stakeholder to drag the other particle

#ub_e = upper_bound_factor*v0
#bounds = ((-ub_e,ub_e))
bounds = ((0,length)) #move this to parameters file

#guidance vector parameters
little_r = 0 
k_far = 2.5 #just changed this from 3.0. Does it affect things? #2.0 was too small. 
big_R = k_far*(2*particle_radius)    
G_path = 1
H_path = 0
G_obstacle = -10
H_obstacle = 1    