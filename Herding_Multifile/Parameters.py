import numpy as np

##Parameters

folder = 'Data' #Folder to save to.
save = False #Set to True if you want to save the data as CSV files. 
filename = 'PositionTarget' #File to save to.

#simulation parameters
#FIXME Now that I use a while loop, it would make more sense to create these arrays dynamically. 
n_times = 3001      #number of simulation time steps
t_end = 3000        #end time in seconds
t = np.linspace(0.0,t_end,n_times) 
delta_t = t[1] - t[0] #Timestep for the controller. 

length = 0.0002   #Domain length of control region in meters. 
n_electrodes = 0   #Number of electrodes #FIXME
boundary_scale = 3 #Number of domain lengths between the simulation boundaries. 
dt = 0.0008       #Timestep for the finite difference simulations.

#target function parameters
#target_type must be one of 'position', 'circle', 'figure8', or 'square'.
#Note that the default initial position for circle is right side, while default for figure8 is center.
#A 'position' target is the only type that can handle more than one passive particle. 
target_type = 'position' 
v0 = 1E-6                               #loiter velocity, m/s
v0_stakeholder = 20*v0 #this is supposed to be the speed the stakeholder moves at
xt = length/2                           #center of the circle
yt = length/2
rd = 3.0E-5                             #desired loiter distance, m
# initial_position = np.array([[xt+2.0*rd,yt+1.0*rd],[xt+rd,yt+2*rd],[xt-rd,yt-rd],[xt+rd,yt-2*rd],[xt+rd,yt],[0,0],[0,0],[0,0],[0,0]]) #initial positions of particles
initial_position = np.zeros((7,2))  #The initial positions of particles (both passive and stakeholders).
for i in range(len(initial_position)):
    initial_position[i] = length/4+np.random.rand(2)*length/2 #for testing
#The first n_stakeholders entries are the initial positions of the stakeholders.
n_particles = len(initial_position)   #The number of particles (both passive and stakeholders).

target_position = np.array([[0,0],[0.0001,0.00012],[0.00012,0.00012],[0.0001,0.00008],[0.00012,0.00008],[0.00008,0.00012],[0.00008,0.00008]])
n_stakeholders = 1                       #number of particles designated as stakeholders. Stakeholders are the first elements of initial_position.

#physical parameters
mu = -2.0E-10/1000                    #diffusiophoretic mobility of particle m^2/M s. Divide by 1000 to turn Liters to meters^3.                           
mu_e_stakeholder = 2*2.0E-8          #electrophoretic mobility of particles in Coulomb/Newton * m/s or m^2/V s
epsilon = 78.4*8.8542E-12            #Dielectric constant of water in C^2/Nm^2 
#diffusion coefficient of solute.
D = 1.38E-9 #urea in water at 25C. 
brownian_motion = 1                            #Brownian motion. Set to 1 for on or 0 for off
T = 298                              #Temperature in K
viscosity = 8.9E-4                   #viscosity of water
particle_radius = 3E-6                        #m, radius of colloidal particle
k_boltzmann = 1.38064852E-23 #Boltzmann's constant.
Dp = k_boltzmann*T/(viscosity*6*np.pi*particle_radius) #The diffusion coefficient of the colloidal particles.

#controller parameters
n_vision = 3#was 5         #prediction horizon
W_1 = 20             #weight for residuals in objective function. This is actually W1^2, the way it is written in the paper.
W_2 = 1
W_3 = 1000
cutoff_percent = 0.2 #higher value will solve faster but be less accurate.
feedback_gain = 0.1  #gain for target function
k_close = 1.1
upper_bound_factor = 50       #multiple of v0_stakeholder to set upper bound at.
#strength = 0.5*2*particle_radius/length #I think this would be the strength needed for the stakeholder to drag the other particle
ube = upper_bound_factor*v0_stakeholder

bounds = ((0,length)) #move this to parameters file

#guidance vector parameters
little_r = 0 
k_far = 2.5 #just changed this from 3.0. Does it affect things? #2.0 was too small. 
big_R = k_far*(2*particle_radius)    
G_path = 1
H_path = 0
G_obstacle = -10
H_obstacle = 1    




#some more setup stuff
decision_times = t.copy() + 0.0000001    #array for decision times. 0.0000001 lets us avoid dividing by zero. 


u_gain = 10
umax = 8*np.pi*D*Dp/np.abs(mu)*u_gain   #characteristic scale for chemical reaction


simulation_length = boundary_scale*length  # Length of simulation domain.
nxy     = 50*boundary_scale+1  # number of grid points in x and y directions
tFac    = 0.5                 # Stability factor for FTCS method. Must be 0.5 or less.
diffusive_timescale     = simulation_length**2/D            # diffusive timescale based on the whole domain.
deltaxy = simulation_length/(nxy-1)         # grid spacing in x and y directions
dtwant = tFac*deltaxy**2/D/4  # desired timestep size
if dt > dtwant:               #make sure to set this so that dt is a nice number smaller than dtwant
    print(dtwant)
    raise
dtratio = int(delta_t/dt)

