import numpy as np
from scipy.optimize import linear_sum_assignment

##Parameters

folder = 'Data' #Folder to save to.
save = False #Set to True if you want to save the data as CSV files. #The save function is not yet implemented
filename = 'PositionTarget' #File to save to.

#simulation parameters
#FIXME Now that I use a while loop, it would make more sense to create these arrays dynamically. 
# =============================================================================
# n_times = 3001      #number of simulation time steps
# t_end = 1000        #end time in seconds
# t = np.linspace(0.0,t_end,n_times) 
# delta_t = t[1] - t[0] #Timestep for the controller. 
# =============================================================================
delta_t = 0.1
t = np.array([0,delta_t])

length = 0.0002   #Domain length of control region in meters. 
n_electrodes = 0   #Number of electrodes #FIXME
boundary_scale = 3 #Number of domain lengths between the simulation boundaries. 
dt = 0.0008       #Timestep for the finite difference simulations.

#target function parameters
#target_type must be one of 'position', 'circle', 'figure8', or 'square'.
#Note that the default initial position for circle is right side, while default for figure8 is center.
#A 'position' target is the only type that can handle more than one passive particle. 
target_type = 'position' 
v0_stakeholder = 4*5*1E-6 #The speed the stakeholder moves at. #was previously 20
xt = length/2                           #center of the circle
yt = length/2

#was 7
initial_position = np.zeros((10,2))  #The initial positions of particles (both passive and stakeholders).
for i in range(len(initial_position)):
    initial_position[i] = length/4+np.random.rand(2)*length/2 #for testing
#The first n_stakeholders entries are the initial positions of the stakeholders.
n_particles = len(initial_position)   #The number of particles (both passive and stakeholders).

center = length/2
ell_target = 0.00004
target_hi = center+ell_target
target_lo = center-ell_target
#target_position_unsorted = np.array([[0,0],[0.0001,0.00012],[0.00012,0.00012],[0.0001,0.00008],[0.00012,0.00008],[0.00008,0.00012],[0.00008,0.00008]])
#target_position_unsorted = np.array([[0,0],[0.0001,0.00014],[0.00014,0.00014],[0.0001,0.00006],[0.00014,0.00006],[0.00006,0.00014],[0.00006,0.00006]])
target_position_unsorted = np.array([[0,0],[center,target_hi],[target_hi,target_hi],[center,target_lo],[target_hi,target_lo],[target_lo,target_hi],[target_lo,target_lo],[target_lo,center],[center,center],[target_hi,center]])
n_stakeholders = 1  #change this to n_herders. They're not stakeholders any more.

#Assign target positions to particles.
distances = np.zeros((n_particles-n_stakeholders,n_particles-n_stakeholders))                     #number of particles designated as stakeholders. Stakeholders are the first elements of initial_position.
for i in range(n_particles-n_stakeholders):
    for j in range(n_particles-n_stakeholders):
        distances[i,j] = np.sqrt((initial_position[i+n_stakeholders,0]-target_position_unsorted[j+n_stakeholders,0])**2 + (initial_position[i+n_stakeholders,1]-target_position_unsorted[j+n_stakeholders,1])**2)
particles, targets = linear_sum_assignment(distances)
target_position = target_position_unsorted.copy()
for i,key in enumerate(targets):
    target_position[i+n_stakeholders] = target_position_unsorted[key+n_stakeholders]

#physical parameters
mu = 2.0E-10/1000                    #diffusiophoretic mobility of particle m^2/M s. Divide by 1000 to turn Liters to meters^3.                           
mu_e_stakeholder = 2*2.0E-8          #electrophoretic mobility of particles in Coulomb/Newton * m/s or m^2/V s
epsilon = 78.4*8.8542E-12            #Dielectric constant of water in C^2/Nm^2 
D = 2.01E-9 #Diffusion coefficient of solute. This number is for oxygen in water. 
brownian_motion = 1                            #Brownian motion. Set to 1 for on or 0 for off
T = 298                              #Temperature in K
viscosity = 8.9E-4                   #viscosity of water
particle_radius = 3E-6                        #m, radius of colloidal particle
herder_radius = 1.0*particle_radius 
k_boltzmann = 1.38064852E-23 #Boltzmann's constant.
Dp = k_boltzmann*T/(viscosity*6*np.pi*particle_radius) #The diffusion coefficient of the colloidal particles.
Dp_herder = k_boltzmann*T/(viscosity*6*np.pi*herder_radius) #not yet implemented

#controller parameters
n_vision = 3#was 5         #prediction horizon
W_1 = 20             #weight for residuals in objective function. This is actually W1^2, the way it is written in the paper.
W_2 = 1
W_3 = 1000
cutoff_percent = 0.2 #higher value will solve faster but be less accurate.
#k_close = 1.1
R_close = 0.1*particle_radius

bounds = ((0,length)) #move this to parameters file

#guidance vector parameters
little_r = 0 
k_far = 1.5#2.5 #just changed this from 3.0. Does it affect things? #2.0 was too small. 
big_R = k_far*(particle_radius + herder_radius)    
G_path = 1
H_path = 0
G_obstacle = -2
H_obstacle = 0.5    

d_precision = 2*particle_radius


#some more setup stuff
decision_times = t.copy() + 1e-7    #array for decision times. 1e-7 lets us avoid dividing by zero. 


reaction_rate = 0.01 #mol/s/m^2, typical of platinum in H2O2.
umax = reaction_rate*2*np.pi*herder_radius
umax = 4e-7


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

#%%
#print(np.sqrt(2*Dp*delta_t))
print(R_close)
vmax = umax*mu/(2*np.pi*D*(R_close+particle_radius+herder_radius))
print(vmax*delta_t)
