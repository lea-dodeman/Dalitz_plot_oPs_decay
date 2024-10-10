import math as m
import numpy as np
from scipy.optimize import fsolve # Non-linear solver
import pandas as pd 
import matplotlib.pyplot as plt

# Conditions :
# E1 + E2 + E3 = MP (total energy conservation)
# theta_12 + theta_13 + theta_23 = 360  (total momentum conservation)
# E3 + E2*cos(theta_23) + E1*cos(theta_23+theta_12) = 0 (momentum conservation in x direction)
# E2*sin(theta_23) + E1*sin(theta_23+theta_12) = 0 (momentum conservation in y direction)

# Constants
MP = 1022  # keV mass of Positronium
EPSILON = 1 # Tolerance for floating-point comparisons
STEP = 0.5  # Step size for angles in degrees
E1_res = 5,4 #keV (Energy resolution for detector 1)
E2_res = 5,7 #keV (Energy resolution for detector 2)
E3_res = 3,8 #keV (Energy resolution for detector 3)


def energy_calc(theta12, theta13, theta23):

    """ Calculate the energies of the three gammas in the decay of Positronium into 3 photons, 
    given the angles between the momenta. """

    # Convert angles to radians
    theta12 = m.radians(theta12)
    theta13 = m.radians(theta13)
    theta23 = m.radians(theta23)

    # Define the system of equations based on special relativity 
    def equations(E):
        E_1, E_2, E_3 = E
        eq1 = E_1 + E_2 + E_3 - MP # Total energy conservation
        eq2 = E_3 + E_2 * m.cos(theta23) + E_1 * m.cos(theta23 + theta12) # Momentum conservation in x direction
        eq3 = E_2 * m.sin(theta23) + E_1 * m.sin(theta23 + theta12) # Momentum conservation in y direction
        return eq1, eq2, eq3
    
    # Initial guess for the energies
    E_initial = [MP / 3, MP / 3, MP / 3]

    # Solve the system of equations
    E_1, E_2, E_3 = fsolve(equations, E_initial)

    return E_1, E_2, E_3

# Function to calculate partial decay rate for positronium decay (formula from the PDG booklet)
def decay_rate_E1(E_1, E_2, E_3, M, matrix_element):
    """ Calculate the partial decay rate for the decay of Positronium into 3 photons
    Comment on the matrix element squared : calculated from QED theory, 
    here supposed to be constant because we consider only 1 particular decay. 
    We will decide to set it to 1, because we normalize the decay rate later."""
    return (1 / (8 * M * (2 * np.pi)**3)) * matrix_element * E_1 * E_2 * E_3

# Create a list of possible configurations of angles
# Must be between 0 and 180 degrees, and must satisfy the condition theta_12 + theta_13 + theta_23 = 360
configs = []
for theta_12 in np.arange(1, 180, STEP):
    for theta_13 in np.arange(1, 180, STEP) :
        theta_23 = 360 - theta_12 - theta_13
        if 0 < theta_23 < 180:
            configs.append([theta_12, theta_13, theta_23])

#Calculate energies and store the results
results = []
for i, config in enumerate(configs):
    theta_12, theta_13, theta_23 = config
    if abs(360 - (theta_12 + theta_13 + theta_23)) > EPSILON : # momentum conservation
        results.append({
            "Configuration": i + 1,
            "theta_12": theta_12,
            "theta_13": theta_13,
            "theta_23": theta_23,
            "E1": "Invalid",
            "E2": "Invalid",
            "E3": "Invalid",
            "Verification": "Invalid"
        })
    else:
        E1, E2, E3 = energy_calc(theta_12, theta_13, theta_23)
        if E1 < 0 or E2 < 0 or E3 < 0 or E1 > 511 or E2 > 511 or E3 > 511: # energy conservation
            results.append({
                "Configuration": i + 1,
                "theta_12": theta_12,
                "theta_13": theta_13,
                "theta_23": theta_23,
                "E1": "Invalid",
                "E2": "Invalid",
                "E3": "Invalid",
                "Verification": "Invalid"
            })
        else:
            E1, E2, E3 = energy_calc(theta_12, theta_13, theta_23)
            m12_squared = (E1 + E2)**2 - (E1**2 + E2**2)
            m13_squared = (E1 + E3)**2 - (E1**2 + E3**2)
            m23_squared = (E2 + E3)**2 - (E2**2 + E3**2)
            dGamma = decay_rate_E1(E1, E2, E3, MP, 1)
            results.append({
                "Configuration": i + 1,
                "dGamma": dGamma,
                "theta_12": theta_12,
                "theta_13": theta_13,
                "theta_23": theta_23,
                "E1": E1,
                "E2": E2,
                "E3": E3,
                "m12_squared": m12_squared,
                'm13_squared': m13_squared,
                "m23_squared": m23_squared,
                "Verification": E1 + E2 + E3
            })

print('len(results) : ', len(results))

df = pd.DataFrame(results)

# Filter out invalid entries
valid_df = df[df['dGamma'] != "Invalid"]

# Calculate the range of decay rates
min_dGamma = valid_df['dGamma'].min()
max_dGamma = valid_df['dGamma'].max()

# Useful prints
#print(f"Decay rate range : {min_dGamma} to {max_dGamma}")
#print('Configuration with highest decay rate : ', valid_df.loc[valid_df['dGamma'].idxmax()])

# Normalize the decay rates
valid_df['dGamma'] = (valid_df['dGamma'] - min_dGamma) / (max_dGamma - min_dGamma)

# Save the DataFrame to a data file
# file_path = 'C:\\Users\\leado\\Documents\\M2_NPAC\\valid_df.pkl'  # Replace with your desired file path
# valid_df.to_pickle(file_path)

# Create a canva for 2 separate subplots
fig, axs = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [1, 1.5]})

# Histogram of the energy E1
sc0 = axs[0].hist(valid_df['E1'], weights=valid_df["dGamma"], bins=100, color='white', edgecolor='purple')
axs[0].set_xlabel('Energy E1 (keV)')
axs[0].set_ylabel('Number of counts')
axs[0].set_title(r'Theoretical count proportion for each photon energy $E_1$')

# Dalitz plot with angles (theta_12, theta_23)
sc1 = axs[1].scatter(valid_df['theta_12'], valid_df['theta_23'], c=valid_df['dGamma'], cmap='magma', s=1)
axs[1].set_xlabel(r'$\theta_{12}$ (degrees)')
axs[1].set_ylabel(r'$\theta_{23}$ (degrees)')
axs[1].set_title(r'Dalitz plot prediction of kinematics for o-Ps decay into 3$\gamma$')
axs[1].set_xlim(-10, 190)
axs[1].set_ylim(-10, 190)
cbar1 = fig.colorbar(sc1, ax=axs[1]) #
cbar1.set_label('Partial Decay Rate dΓ')

# Find the configuration with the highest decay rate
highest_decay_rate_config = valid_df.loc[valid_df['dGamma'].idxmax()]
print('Configuration with highest decay rate : ', highest_decay_rate_config)

# Add annotation for the configuration with the highest decay rate
highest_decay_text = (f"Highest probability for :\n"
                       f"$\\theta_{12}$ = $\\theta_{23}$ = 120°")
axs[1].text(0.05, 0.05, highest_decay_text, transform=axs[1].transAxes, fontsize=12, 
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.show()
