import matplotlib.pyplot as plt
import numpy as np

# Input parameters (in centimeters)
length_cm = 5                # Length of the cylindrical system in cm
diameter_cathode_cm = 0.12   # Diameter of the cathode in cm
diameter_anode_cm = 2        # Diameter of the anode in cm
diameter_grid_cm = 0.4       # Diameter of the grid in cm
grid_spacing_cm = 0.25       # Spacing between grid wires in cm
diameter_grid_wire_cm = 0.01 # Diameter of the grid wire in cm

# Emission ratio
emr = 0.03  # A / cm^2

# Calculation of the cathode area S
r_cathode = diameter_cathode_cm / 2  # Radius of the cathode in centimeters
S = 2 * np.pi * r_cathode * length_cm  # Cathode area in cm^2

# Calculation of the saturation current from the area and emission ratio
I_saturation = emr * S

# Calculation of radii and ratio r/a
a      = diameter_cathode_cm / 2    # Radius of the cathode in cm
r      = diameter_anode_cm / 2      # Radius of the anode in cm
r_grid = diameter_grid_cm / 2       # Radius of the grid in centimeters
rg     = diameter_grid_wire_cm / 2  # Radius of the grid wire in centimeters
r_a    = r / a                      # Ratio r/a
l      = length_cm                  # Length of electrodes (height of the cylinder)

# Calculation of the amplification factor μ
xa    = rho_p = r                # Radius of the anode in cm
xg    = rho_g = r_grid           # Radius of the grid in centimeters
p     = a     = grid_spacing_cm  # Spacing between grid wires (pitch) in cm
#rg                              # Radius of the grid wire in cm

# Kusunose simplified spiral length
Lg1 = 2 * np.pi * xg / p

# Kusunose exact spiral length
Lg2 = np.sqrt(1 + (Lg1)**2)

# Selection of Lg
Lg = Lg1

# Kusunose C factor
C = 1

# Reich (3 - 5)
mu  = - (2 * np.pi * (rho_g / rho_p) * (rho_p - rho_g)) / (a * np.log(a / (2 * np.pi * rg)))

# Kusunose (k - 2)
mu2 = - Lg * (np.log(xa / xg) / np.log(p / (2 * np.pi * rg)))

# Kusunose (k - 3)
mu3 = - (Lg * np.log(xa / xg) - np.log(np.cosh(2 * np.pi * rg / p))) / (np.log(np.cosh(2 * np.pi * rg / p)) - np.log(np.sinh(2 * np.pi * rg / p)))

# Kusunose (k - 4)
mu4 = - C * Lg * (np.log(xa / xg) / np.log(1 / np.tanh(2 * np.pi * rg / p)))

# Kusunose final eq. 30 (k - 5)
Pia = 2 * np.pi* rg / p
mu5 = - C * Lg * (np.log(xa / xg) / np.log(1 / np.tanh(Pia)))

print(f"Amplification factor μ (Reich): {mu:.2f}, (k-2): {mu2:.2f}, (k-3): {mu3:.2f}, (k-4): {mu4:.2f} (k-5): {mu5:.2f}")

# Selection of mu for further calculations
mu_abs = np.abs(mu5)

# Calculation of beta for a given r/a using the specified formula
log_ra = np.log(r_a)
beta = (log_ra - 
        (2/5) * (log_ra)**2 + 
        (11/120) * (log_ra)**3 - 
        (47/3300) * (log_ra)**4)

# Range of voltage values V in volts
V_values = np.linspace(1, 400, 50)

# Range of grid voltage values Ug
Ug_values = np.arange(-20, 20, 10)  # From -20V to 20V in steps of 10V

# Plotting the graph for each Ug value
plt.figure(figsize=(10, 6))

for Ug in Ug_values:
    # Calculation of current (Child-Langmuir Law) for a cylindrical system using the formula in SI units (be careful on e.s.u. units) - adjusted for a triode
    #I_values = (14.68e-6 * l /(r * beta**2)) * ((np.maximum(0, V_values + Ug * mu_abs))**1.5)
    
    # Taking mu into account according to Kusunose
    #I_values = (14.68e-6 * l /(r * beta**2)) * ((np.maximum(0, V_values + Ug * mu_abs))**1.5) / ((1 + mu_abs)**1.5)

    # Full calculation according to Kusunose
    ea = V_values
    eg = Ug
    k = mu_abs
    A = 2 * np.pi * xa * l
    G = 2.33e-6 * A / (xa * xg)
    eg_ = (ea + k * eg) / (1 + k)
    I_values = G * (np.maximum(0, eg_))**1.5

    # Attempt to account for saturation
    I_space_sat_values = I_saturation * np.tanh(I_values / I_saturation)

    I_values = I_values * 1e3  # Conversion to mA
    I_space_sat_values = I_space_sat_values * 1e3  # Conversion to mA

    # Plotting the curve for the given Ug
    plt.plot(V_values, I_space_sat_values, label=f'Ug = {Ug} V')

I_sat_values = I_saturation * np.ones_like(V_values)
I_sat_values = I_sat_values * 1e3  # Conversion to mA
plt.plot(V_values, I_sat_values, label=f'Saturation = {I_saturation * 1e3:.1f} mA', linestyle=':', color='tab:green')

# Configuring the graph
plt.xlabel('Anode voltage (Ua) [V]')
plt.ylabel('Anode current (Ia) [mA]')
plt.title('Comparison of currents I at voltage V for different Ug values')
plt.grid(True)
plt.legend()

# Displaying the graph
plt.show()
