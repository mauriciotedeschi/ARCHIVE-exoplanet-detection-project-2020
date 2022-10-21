import math
import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.674e-11 # Gravitational constant, m^3 kg^-1 s^-2
MASS_SUN = 1.989e30 # kilograms
MASS_EARTH = 5.972e24 # kilograms
METER_TO_AU_RATIO = 1.495978707e11
SECOND_TO_DAY_RATIO = 8.6400e4
SPEED_OF_LIGHT = 299792458.0 # m/sec

def convert_polar_to_rect(r, theta, phi, degrees=False):
    if degrees:
        theta = theta * math.pi / 180.0
        phi = phi * math.pi / 180.0
        
    x = r * math.cos(phi) * math.cos(theta)
    y = r * math.cos(phi) * math.sin(theta)
    z = r * math.sin(phi)
    return [x, y, z]

# Calculating doppler shift from velocity
def get_ds_from_vel(v):
    beta = v / SPEED_OF_LIGHT
    return math.sqrt((1+beta)/(1-beta)) - 1
    
# Note: centered at 0
def get_vel_from_ds(z):
    k = (1 + z)**2
    return SPEED_OF_LIGHT * (k - 1) / (k + 1)

# Let x,y dimensions be coplanar to the computer screen, and z dimension projected toward the observer.

def euler_solution(mass_star, pos_star, vel_star, mass_planet, pos_planet, vel_planet, n_steps=100, total_time=86400.0):
    # Conversion of input vectors into NumPy arrays
    initial_planet_position = np.array(pos_planet)
    initial_planet_velocity = np.array(vel_planet)
    initial_star_position = np.array(pos_star)
    initial_star_velocity = np.array(vel_star)
    
    # Initialize lists
    time = np.zeros(n_steps+1)
    pos = np.zeros((2, n_steps+1, 3))
    vel = np.zeros((2, n_steps+1, 3))
    
    pos[0, 0] = initial_star_position
    pos[1, 0] = initial_planet_position
    vel[0, 0] = initial_star_velocity
    vel[1, 0] = initial_planet_velocity
    
    delta_t = total_time / n_steps # Seconds
       
    # Start integration
    for i in range(n_steps):
        r1 = pos[0, i] # Current star position vector
        r2 = pos[1, i] # Current planet posiiton vector
        dist = math.sqrt(((r1 - r2)**2).sum()) # Distance between star and planet
        F_g = G * mass_star * mass_planet / dist**2
        
        unit_vector_star = (pos[1, i] - pos[0, i]) / dist
        unit_vector_planet = -unit_vector_star
        acc_p = F_g / mass_planet * unit_vector_planet # Acceleration vector for planet
        acc_s = F_g / mass_star * unit_vector_star # Acceleration vector for star
        
        vel[0, i + 1] = vel[0, i] + acc_s * delta_t
        vel[1, i + 1] = vel[1, i] + acc_p * delta_t
        pos[0, i + 1] = pos[0, i] + vel[0, i + 1] * delta_t
        pos[1, i + 1] = pos[1, i] + vel[1, i + 1] * delta_t
        time[i + 1] = time[i] + delta_t
    return [time, pos, vel]

def plot_euler_solution(t, r, v, figsize): # Outputs a figure object that can be saved or displayed
    fig, ax = plt.subplots(1,2,figsize=figsize)
    ax[0].plot(r[1,:,0],r[1,:,2], 'b') # Plot orbits of planet and star in first subplot
    ax[0].plot(r[0,:,0],r[0,:,2], 'y')
    ax[0].set(title='Positions of Planet and Star', xlabel='x (m)', ylabel ='z (m)')
    #ax[0].set_xlim(-5e11, 5e11)
    #ax[0].set_ylim(-10e11, 10e11)
    
    ax[1].plot(t, v[0,:,2]) # Plot star velocity parallel to "observer" vs time
                            # Let y-component be parallel to line of sight
    ax[1].set(title='Stellar Radial Velocity', xlabel='Time (seconds)', ylabel='Velocity to observer (ms^-1)')
    
    return fig, ax

def save_euler_solution_to_file(t, r, v, figsize=(6,4)):
    my_plot, ax = plot_euler_solution(t, r, v, figsize)
    i = 0
    try:
        i = int(open('nbody_filetracker.txt','r').readline(-1))
        open('nbody_filetracker.txt','w').write(str(i + 1))        
    except:
        f = open('nbody_filetracker.txt','w')
        f.write('0')
        f.close()
    file_out = 'plot_' + str(i) + '.png'
    my_plot.savefig(file_out)
    
def save_plot_to_file(fig, ax):
    i = 0
    try:
        i = int(open('nbody_filetracker.txt','r').readline(-1))
        open('nbody_filetracker.txt','w').write(str(i + 1))        
    except:
        f = open('nbody_filetracker.txt','w')
        f.write('0')
        f.close()
    file_out = 'plot_' + str(i) + '.png'
    fig.savefig(file_out)
    
def main1():
    print('Testing nbody simulation')
    
    t, r, v = euler_solution(MASS_SUN, [0.0,0.0,0.0], [-2978.8,0.0,0.0], 
                                0.1 * MASS_SUN, [0.0, 0.0, METER_TO_AU_RATIO], [29788.0, 0.0,0.0], 
                                total_time=100*SECOND_TO_DAY_RATIO, n_steps=300)
    print(t[-1])
    save_euler_solution_to_file(t, r, v, figsize=(12,12))
    print('Done')
    
def main2():
    sol = euler_solution(MASS_SUN, [0.0,0.0,0.0], [-2978.8,0.0,0.0], 
                                0.1 * MASS_SUN, [0.0, METER_TO_AU_RATIO,0.0], [29788.0, 0.0,0.0], 
                                total_time=100*SECOND_TO_DAY_RATIO, n_steps=300)
    print('Star velocity', sol[1][0,:,1])
    
if __name__ == "__main__": # Test to see that everything works
    main1()




