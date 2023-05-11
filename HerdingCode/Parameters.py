import numpy as np
from scipy.optimize import linear_sum_assignment

##Parameters

folder = 'Data' #Folder to save to.
save = True #Set to True if you want to save the data as CSV files. #The save function is not yet implemented
filename = 'PositionTarget' #File to save to.

#simulation parameters
delta_t = 0.1
t = np.array([0,delta_t])

length = 0.0002   #Domain length of control region in meters. 
boundary_scale = 3 #Number of domain lengths between the simulation boundaries. 
dt = 0.0008       #Timestep for the finite difference simulations.


v0_stakeholder = 5E-6 #The speed the stakeholder moves at. #was previously 20

initial_position = np.zeros((9,2))  #The initial positions of particles (both passive and stakeholders).
for i in range(len(initial_position)):
    initial_position[i] = length/4+np.random.rand(2)*length/2 #for testing
#The first n_herders entries are the initial positions of the stakeholders.
n_particles = len(initial_position)    #The number of particles (both passive and stakeholders).
#If I were to rewrite this, I would have stakeholders and followers in different arrays. The way I have it here is too confusing.

center = length/2
ell_target = 0.00005 #spacing of targets
target_position_unsorted = np.zeros((n_particles,2))
n_herders = 1  

circle_radius = ell_target*(n_particles-1)/(2*np.pi)
for i in range(1,n_particles): #space the targets evenly around the circle
    angle = 2*np.pi/(n_particles-1)*(i-1)
    target_position_unsorted[i] = [center + circle_radius*np.cos(angle),center + circle_radius*np.sin(angle)]


# =============================================================================
# ell_target = 0.00005
# target_hi = center+ell_target
# target_lo = center-ell_target
# target_hi2 = center+ell_target/2
# target_lo2 = center-ell_target/2
# 
# target_position_unsorted = np.array([[0,0],[center,target_hi2],[target_hi,target_hi2],[center,target_lo2],[target_hi,target_lo2],[target_lo,target_hi2],[target_lo,target_lo2]])#,[target_lo,center],[center,center],[target_hi,center]])
# =============================================================================



#Assign target positions to particles.
distances = np.zeros((n_particles-n_herders,n_particles-n_herders))                     #number of particles designated as stakeholders. Stakeholders are the first elements of initial_position.
for i in range(n_particles-n_herders):
    for j in range(n_particles-n_herders):
        distances[i,j] = np.sqrt((initial_position[i+n_herders,0]-target_position_unsorted[j+n_herders,0])**2 + (initial_position[i+n_herders,1]-target_position_unsorted[j+n_herders,1])**2)
particles, targets = linear_sum_assignment(distances)
target_position = target_position_unsorted.copy()
for i,key in enumerate(targets):
    target_position[i+n_herders] = target_position_unsorted[key+n_herders]

#physical parameters
mu = 2.0E-10/1000                    #diffusiophoretic mobility of particle m^2/M s. Divide by 1000 to turn Liters to meters^3.                           
epsilon = 78.4*8.8542E-12            #Dielectric constant of water in C^2/Nm^2 
D = 2.01E-9 #Diffusion coefficient of solute. This number is for oxygen in water. 
brownian_motion = 1                            #Brownian motion. Set to 1 for on or 0 for off
T = 298                              #Temperature in K
viscosity = 8.9E-4                   #viscosity of water
particle_radius = 3E-6                        #m, radius of colloidal particle
herder_radius = 3e-6
k_boltzmann = 1.38064852E-23 #Boltzmann's constant.
Dp = k_boltzmann*T/(viscosity*6*np.pi*particle_radius) #The diffusion coefficient of the follower particles.
Dp_herder = k_boltzmann*T/(viscosity*6*np.pi*herder_radius) #The diffusion coefficient of the herders.

#controller parameters
R_close = 0.2*particle_radius
d_prec = 1e-6 #herding precision
d_tol = 9e-6 #herding tolerance, for when to stop algorithm. need to find an expression for this, it is incomplete.

#guidance vector parameters
little_r = 0 
k_far = 1.6
big_R = k_far*(particle_radius + herder_radius)    
G_path = 1
H_path = 0
G_obstacle = -2
H_obstacle = 0.5    



#some more setup stuff
decision_times = t.copy() + 1e-7    #array for decision times. 1e-7 lets us avoid dividing by zero. 


reaction_rate = 0.008 #mol/s/m^2, 0.01 is typical of platinum in H2O2. So this is 1/3 of the max physically reasonable.
umax = reaction_rate*2*np.pi*herder_radius


simulation_length = boundary_scale*length  # Length of simulation domain.
nxy     = 50*boundary_scale+1  # number of grid points in x and y directions
tFac    = 0.5                 # Stability factor for FTCS method. Must be 0.5 or less.
diffusive_timescale     = simulation_length**2/D            # diffusive timescale based on the whole domain.
deltaxy = simulation_length/(nxy-1)         # grid spacing in x and y directions
dtwant = tFac*deltaxy**2/D/4  # desired timestep size
if dt > dtwant:               #make sure to set this so that dt is a nice number smaller than dtwant
    print(dtwant)
    raise #clarify this
dtratio = int(delta_t/dt)


#run checks on dprec
dbrown1 = np.sqrt(4*Dp*delta_t)
if d_prec < dbrown1:
    print("Brownian precision check failed")

#run check on dtol
kdiff = umax*mu/(2*np.pi*D)
R_fh = particle_radius + herder_radius + R_close
therder = ell_target/v0_stakeholder
tchased =  d_tol*R_fh/kdiff
#ttotal = (n_particles-n_herders)*(tchased+therder) #this n_particles is not the same number used in the paper because this contains stakeholders too.
#dbrown2 = np.sqrt(4*Dp*ttotal)

dbrown2 = 2*Dp*(n_particles-n_herders)*R_fh/kdiff + 2*np.sqrt((Dp*(n_particles-n_herders)*R_fh/kdiff)**2 + Dp*(n_particles-n_herders)*therder)

#dbrown2 = 4*Dp*(n_particles-1)*R_fh/kdiff #this n_particles is not the same number used in the paper because this contains stakeholders too.
if d_tol < dbrown2:
    print("Brownian tolerance check failed, dtol=",d_tol,"dbrown=",dbrown2)


#%%
#vmax = umax*mu/(2*np.pi*D*(R_close+particle_radius+herder_radius))
#print(vmax*delta_t)
#print(kdiff)
#print(umax)
print(d_tol,dbrown2)