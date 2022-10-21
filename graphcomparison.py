import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import rvdatareader
import nbodysolver as nbody

# Numbers we can consider as constants
star_masses = rvdatareader.read_star_mass_data('star-masses.csv')

def get_SST(x_recorded, y_recorded, x_simulated, y_simulated):
    func = interpolate.interp1d(x_simulated, y_simulated, kind = 'cubic')
    # sum of square residuals
    tot = 0
    for i in range(len(x_recorded)):
        y_hat = func(x_recorded[i])
        tot += (y_recorded[i] - y_hat)**2
    return tot
    
# Should give us the time and velocity to start from in Nbody integration.
def get_max_vel_point(t_rec, v_rec):
    t0, max_vel = 0, 0
    for i in range(len(t_rec)):
        if abs(v_rec[i]) > abs(max_vel):
            t0 = t_rec[i]
            max_vel = v_rec[i]
            
    return [t0, max_vel]


def calc_exoplanet_data(filename, star_mass_index, # index: denotes which star we are working with.
                        mass_planet_min_guess, 
                        mass_planet_max_guess, 
                        orbital_radius_min_guess,
                        orbital_radius_max_guess,
                        num_trials=64,
                        num_steps_for_nbody=300):
                        
    # Get star mass from list we compiled above.
    mass_star = star_masses[star_mass_index] * nbody.MASS_SUN
    # Calculate time interval needed for n-body sim
    t_recorded, ds_recorded = rvdatareader.read_doppler_shift_data(filename)
    total_time = t_recorded[-1] - t_recorded[0]
    # We will translate recorded dopplershift into radial velocity
    v_recorded = nbody.get_vel_from_ds(np.array(ds_recorded))
    #v_shift = (max(v_recorded) + min(v_recorded)) / 2.0
    #v_recorded -= v_shift
    
    # We will find the time at which graph hits maximum velocity
    t0, v_star = get_max_vel_point(t_recorded, v_recorded)
    
    best_t_sim, best_v_sim = [], []
    best_r_sim = []
    
    print('Star mass:', mass_star)
    print('Velocity of star:', v_star)
    print('Guessing')
    
    min_sst = np.inf # infinity
    
    for step in range(num_trials):
        
        # Find the geometric means of guess bounds of both variables
        mass_planet = np.sqrt(mass_planet_min_guess * mass_planet_max_guess)
        radius = np.sqrt(orbital_radius_min_guess * orbital_radius_max_guess)
        
        # Construct the solution and get SST
        r1 = np.sqrt(radius * orbital_radius_min_guess)
        r2 = np.sqrt(radius * orbital_radius_max_guess)
        mp1 = np.sqrt(mass_planet * mass_planet_min_guess)
        mp2 = np.sqrt(mass_planet * mass_planet_max_guess)
        
        best_m_p, best_r = 0, 0
        
        for r in [r1, r2]:
            for m in [mp1, mp2]:
                r_planet = mass_star * r / (m + mass_star) # Distance of planet from barycenter
                r_star = m * r / (m + mass_star) # Distance of star from barycenter
                v_planet = np.sqrt(nbody.G * mass_star * r_planet) / r
                # Star starts going up in z dimension, planet down in z dimension
               
                r_star_vector = np.array([-r_star, 0.0, 0.0])
                v_star_vector = np.array([0.0, 0.0, v_star])
                r_planet_vector = np.array([r_planet, 0.0, 0.0])
                v_planet_vector = np.array([0.0, 0.0, -v_planet])
                
                soln = nbody.euler_solution(mass_star, r_star_vector, v_star_vector, m, r_planet_vector, v_planet_vector, num_steps_for_nbody, total_time)
                # Recall that the remaining time from 0 to the starting point of former soln is t0
                soln_reverse = nbody.euler_solution(mass_star, r_star_vector, -v_star_vector, m, r_planet_vector, -v_planet_vector, num_steps_for_nbody, t0)
                
                # We can concatenate the two arrays
                # Just remove the first time step in the forward sim so we don't overlap
                t_sim = np.concatenate((soln_reverse[0], soln[0][1:] + t0))
                # For the velocity (its component) we will have to flip the array so that it maps backward correctly
                v_sim = np.concatenate((-np.flip(soln_reverse[2], 1), soln[2][:,1:]), 1)
                
                    
                #print(soln_reverse[0])
                #print(soln[0])
                #print(v_sim[0,:,2])
                
                sst = get_SST(t_recorded, v_recorded, t_sim, v_sim[0,:,2])
                if sst < min_sst:
                    min_sst = sst
                    best_m_p = m
                    best_r = r
                    best_t_sim = t_sim
                    best_v_sim = v_sim
                    best_r_sim = soln_reverse[1]
                    #best_r_sim = np.concatenate((np.flip(soln_reverse[1], 1), soln[1][:,1:]), 1)
                    #best_r_sim = soln[1]
        
        print((step+1), end=' ')
        # Choosing the bounds for the next iteration
        if best_m_p == mp1:
            mass_planet_max_guess = mass_planet
            print('Lower mass', end=' ')
        elif best_m_p == mp2:
            mass_planet_min_guess = mass_planet
            print('Upper mass', end=' ')
        else:
            # Narrow down the bounds in both directions
            mass_planet_min_guess = mp1
            mass_planet_max_guess = mp2
        if best_r == r1:
            orbital_radius_max_guess = radius
            print('Lower radius', end='\t')
        elif best_r == r2:
            orbital_radius_min_guess = radius
            print('Upper radius', end='\t')
        else:
            # Narrow down the bounds in both directions
            orbital_radius_min_guess = r1
            orbital_radius_max_guess = r2
            
        print(min_sst)
    # End loop
    
    return mass_planet, radius, best_t_sim, best_r_sim, best_v_sim
                
                

if __name__ == '__main__':
    # Short test.
    m_p, radius, t, r, v = calc_exoplanet_data('star0.csv', 0, 0.01 * nbody.MASS_EARTH, 10000 * nbody.MASS_EARTH, 
                                    0.01 * nbody.METER_TO_AU_RATIO, 10000 * nbody.METER_TO_AU_RATIO,
                                    num_trials=64, num_steps_for_nbody=300)
    
    rec_time, rec_dop = rvdatareader.read_doppler_shift_data('star0.csv')
    rec_vel = nbody.get_vel_from_ds(np.array(rec_dop))
    print(m_p, radius, sep='\t\t')
    fig, ax = nbody.plot_euler_solution(t, r, v, (16,9))
    ax[1].plot(rec_time, rec_vel, 'o')
    nbody.save_plot_to_file(fig, ax)
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    