import numpy as np

##Parameters
folder = 'data'
save = False #save the data files or not
filename = 'ParamsFigure8'

#simulation parameters
n_times = 601      #number of simulation time steps
t_end = 300        #end time in seconds
length = 0.0002   #domain length in meters. 
n_electrodes = 0   #number of electrodes
boundary_scale = 3 #number of domain lengths extended for finite difference. 3 is one domain and another on each side
dt = 0.0008       #timestep for physics

#target function parameters
target_type = 'figure8' #note that the default start position for circle is right side, while default for figure8 is center.
v0 = 1E-6                               #loiter velocity, m/s
xt = length/2                           #center of the circle
yt = length/2
rd = 3.0E-5                             #desired loiter distance, m
initial_position = np.array([[xt+rd,yt+rd],[xt,yt]]) #initial positions of particles
n_stakeholders = 1                       #number of particles designated as stakeholders. Stakeholders are the first elements of initial_position.

#physical parameters
mu = 2.0E-10/1000                    #diffusiophoretic mobility of particle m^2/M s. Divide by 1000 to turn Liters to meters^3.
mu_e = 0                              #electrophoretic mobility of particles in Coulomb/Newton * m/s or m^2/V s 
mu_stakeholder = 0
mu_e_stakeholder = 2.0E-8
epsilon = 78.4*8.8542E-12            #Dielectric constant of water in C^2/Nm^2 
D = 2.30E-9                          #diffusion coefficient of solute. This number is the self-diffusion coeff of water at 25 C.
Brown = 1                            #Brownian motion. Set to 1 for on or 0 for off
T = 298                              #Temperature in K
kb = 1.38064852E-23                  #Boltzmann's constant
viscosity = 8.9E-4                   #viscosity of water
particle_radius = 3E-6                        #m, radius of colloidal particle

#controller parameters
n_vision = 10         #prediction horizon
W_1 = 20             #weight for residuals in objective function
cutoff_percent = 0.2 #higher value will solve faster but be less accurate.
feedback_gain = 0.1  #gain for target function
ub_e_factor = 10       #multiple of uemax to set upper bound at. 1 is normal for a single particle.
strength = 2.0*2*particle_radius/length #I think this would be the strength needed for the stakeholder to drag the other particle

ub_e = ub_e_factor*v0
bnds = ((-ub_e,ub_e))