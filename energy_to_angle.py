import numpy as np
import scipy.optimize as opt
import math as m
import matplotlib.pyplot as plt

# Constants
mp = 1022  # Mass of ortho-Positronium in keV
d = 55 # Diameter of the detectors in mm (assumed to be the same for all detectors)
u_theta = 1 # Uncertainty in the angle measurement in degrees 
R = 150 # Distance between the source and the detectors in mm

def angular_resolution(theta):

    #  (pi*R/4) * (1 - arctan(dtheta/R))
    return (m.pi * R / 4) * (1 - m.atan(theta/ R))

def energy_calc(theta_a, theta_b, theta_c):
    """ Calculate the energies of the three photons in the decay of Positronium into 3 photons,
    given the angles between the photons and taking into account energy and momentum conservation. """

    theta_a = m.radians(theta_a)
    theta_b = m.radians(theta_b)
    theta_c = m.radians(theta_c)

    # Define the system of equations based on energy and momentum conservation
    def equations(E):
        E_1, E_2, E_3 = E
        eq1 = E_1 + E_2 + E_3 - mp
        eq2 = E_3 + E_2 * m.cos(theta_c) + E_1 * m.cos(theta_c + theta_a)  # Momentum conservation in x direction
        eq3 = E_2 * m.sin(theta_c) + E_1 * m.sin(theta_c + theta_a)  # Momentum conservation in y direction
        return eq1, eq2, eq3
    
    # Initial guess for the energies
    E_initial = [mp / 3, mp / 3, mp / 3]

    # Solve the system of equations
    E_1, E_2, E_3 = opt.fsolve(equations, E_initial)

    return E_1, E_2, E_3

def angles_calc(E1, E2, E3):
    """ Calculate the angles between the three photons in the decay of Positronium into 3 photons,
    given the energies of the photons and taking into account energy and momentum conservation. """
    
    # Define the system of equations based on energy and momentum conservation
    def equations(theta):
        theta_ab, theta_bc, theta_ac = theta
        eq1 = theta_ab + theta_bc + theta_ac - 360  # Total energy conservation
        eq2 = E3 * m.cos(m.radians(theta_ab + theta_bc)) + E2 * m.cos(m.radians(theta_ab)) + E1  # Momentum conservation in x direction
        eq3 = E3 * m.sin(m.radians(theta_ab + theta_bc)) + E2 * m.sin(m.radians(theta_ab))  # Momentum conservation in y direction
        return eq1, eq2, eq3

    # Initial guess for the angles (in degrees)
    theta_initial = [120, 120, 120]

    # Solve the system of equations
    theta_ab, theta_bc, theta_ac = opt.fsolve(equations, theta_initial)

    # Ensure angles are within the range [0, 180] degrees
    theta_ab = theta_ab % 360
    theta_bc = theta_bc % 360
    theta_ac = theta_ac % 360

    if theta_ab > 180:
        theta_ab = 360 - theta_ab
    if theta_bc > 180:
        theta_bc = 360 - theta_bc
    if theta_ac > 180:
        theta_ac = 360 - theta_ac

    return theta_ab, theta_bc, theta_ac

if __name__ == "__main__":

    # Example usage
    E1, E2, E3 = 273, 251.4, 505.1  # Example energies in keV
    theta_ab, theta_bc, theta_ac = angles_calc(E1, E2, E3)
    print(f"Theta_ab: {theta_ab:.2f} degrees, Theta_bc: {theta_bc:.2f} degrees, Theta_ac: {theta_ac:.2f} degrees")

    print(theta_ab + theta_bc + theta_ac) # Total must be 360 degrees



    # Create a list of possible configurations of angles
    configs = []
    for theta_12 in np.arange(1, 180, 0.5):
        for theta_13 in np.arange(1, 180, 0.5) :
            theta_23 = 360 - theta_12 - theta_13
            if 0 < theta_23 < 180:
                configs.append([theta_12, theta_13, theta_23])

    # Filter configurations to keep only those where at least two angles are between 100° and 140°
    filtered_configs = [config for config in configs if sum(110 <= angle <= 130 for angle in config) >= 2]

    print(f"Number of possible configurations: {len(configs)}")
    print(f"Number of filtered configurations: {len(filtered_configs)}")


    # Calculate the energies and store the results
    result_test = []
    fluc = 10
    configuration_test = [120, 120, 120]


    # Calculate the energies for the given configuration
    E1, E2, E3 = energy_calc(*configuration_test)
    result_test = [[E1, E2, E3]]

    # Apply the fluctuation while keeping the sum of angles 360°
    for delta in [-fluc, fluc]:
        fluctuated_config = [angle + delta for angle in configuration_test]
        # Adjust the last angle to ensure the sum is 360°
        #fluctuated_config[-1] = 360 - sum(fluctuated_config[:-1])
        #print(f"Fluctuated configuration: {fluctuated_config}")  # Debug print
        E1, E2, E3 = energy_calc(*fluctuated_config)
        print(f"Calculated energies: E1={E1}, E2={E2}, E3={E3}")  # Debug print
        print(f"Sum of energies: {E1 + E2 + E3}")  # Debug print
        result_test.append([E1, E2, E3])


    # Convert results to a NumPy array
    result_test = np.array(result_test)

    # Calculate the total possible fluctuation in energy
    E1_min, E1_max = np.min(result_test[:, 0]), np.max(result_test[:, 0])
    E2_min, E2_max = np.min(result_test[:, 1]), np.max(result_test[:, 1])
    E3_min, E3_max = np.min(result_test[:, 2]), np.max(result_test[:, 2])

    E1_range = E1_max - E1_min
    E2_range = E2_max - E2_min
    E3_range = E3_max - E3_min

    #print(f"E1 range: {E1_range} keV")
    #print(f"E2 range: {E2_range} keV")
    #print(f"E3 range: {E3_range} keV")

