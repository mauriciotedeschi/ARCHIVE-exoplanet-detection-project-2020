import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import csv

# Returns an array containing the data points of radial velocity vs time
def read_doppler_shift_data(filename):
    time, doppler_shift = [], []
    with open(filename) as file:
        reader = csv.reader(file, delimiter=',')
        line_count = 0
        for row in reader:
            if line_count > 0:
                time.append(float(row[0]))
                doppler_shift.append(float(row[1]))
            line_count += 1
    return time, doppler_shift
    

# This will take in CSV's with data for guess bounds for the guess and check algorithm
def read_guessandcheck_data(filename):
    pass
    """
        Structure:
        will do later
    """
    
# 1D array of star masses denoted by indices
# NOTE: units are Solar masses
def read_star_mass_data(filename):
    masses = []
    with open(filename) as file:
        reader = csv.reader(file, delimiter=',')
        line_count = 0
        for row in reader:
            if line_count > 0:
                masses.append(float(row[1]))    # We only need to take in the second column
            line_count += 1
    return masses
    
if __name__ == "__main__":
    star_masses = read_star_mass_data('star-masses.csv')
    for m in star_masses:
        print(m)