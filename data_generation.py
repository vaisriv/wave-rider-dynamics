import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters
mach_numbers = np.linspace(0.1, 10, 100)  # Mach numbers from 0.1 to 10
angles_of_attack = np.linspace(-5, 15, 21)  # Angles of attack from -5 to 15 degrees

# Aerodynamic coefficients
C_l0 = 0.2    # Lift coefficient at zero AoA
C_l_alpha = 0.1  # Lift curve slope per degree
C_d0 = 0.02   # Zero-lift drag coefficient
k = 0.04      # Induced drag factor

# Generate data
data = []
for mach in mach_numbers:
    for alpha in angles_of_attack:
        alpha_rad = np.radians(alpha)
        if mach < 0.9:
            # Subsonic regime, Prandtl-Glauert correction
            beta = np.sqrt(1 - mach**2)
            C_l_incompressible = C_l0 + C_l_alpha * alpha
            C_l = C_l_incompressible / beta
            C_d = C_d0 + k * C_l**2
        elif 0.9 <= mach < 1.2:
            # Transonic regime, Karman-Tsien correction
            C_l_incompressible = C_l0 + C_l_alpha * alpha
            numerator = C_l_incompressible
            denominator = np.sqrt(1 - mach**2) + (mach**2 * (1 + (mach**2)) * (alpha_rad)**2) / (2 * (1 - mach**2))
            C_l = numerator / denominator
            # Estimate drag divergence
            drag_divergence = 0.1 * np.exp(20 * (mach - 1))
            C_d = C_d0 + k * C_l**2 + drag_divergence
        else:
            # Supersonic regime, linearized supersonic theory
            beta = np.sqrt(mach**2 - 1)
            C_l = 4 * alpha_rad / beta
            # Wave drag
            C_d_wave = (4 * alpha_rad)**2 / (np.pi * beta**2)
            C_d = C_d0 + C_d_wave
        
        # Append to data list
        data.append([mach, alpha, C_l, C_d])

# Create DataFrame
df = pd.DataFrame(data, columns=['Mach', 'AoA', 'C_l', 'C_d'])

# Save to CSV
df.to_csv('aerodynamic_data.csv', index=False)

# Plotting C_l vs AoA
plt.figure(figsize=(10, 6))
for mach in mach_numbers[::10]:  # Plot every 10th Mach number to reduce clutter
    subset = df[df['Mach'] == mach]
    plt.plot(subset['AoA'], subset['C_l'], label=f'Mach {mach:.1f}')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Lift Coefficient (C_l)')
plt.title('Lift Coefficient vs. Angle of Attack for Various Mach Numbers')
plt.legend()
plt.grid(True)
plt.show()

# Plotting C_d vs AoA
plt.figure(figsize=(10, 6))
for mach in mach_numbers[::10]:
    subset = df[df['Mach'] == mach]
    plt.plot(subset['AoA'], subset['C_d'], label=f'Mach {mach:.1f}')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Drag Coefficient (C_d)')
plt.title('Drag Coefficient vs. Angle of Attack for Various Mach Numbers')
plt.legend()
plt.grid(True)
plt.show()

# Plotting C_l / C_d vs AoA
plt.figure(figsize=(10, 6))
for mach in mach_numbers[::10]:
    subset = df[df['Mach'] == mach]
    plt.plot(subset['AoA'], subset['C_l'] / subset['C_d'], label=f'Mach {mach:.1f}')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Lift-to-Drag Ratio (C_l / C_d)')
plt.title('Lift-to-Drag Ratio vs. Angle of Attack for Various Mach Numbers')
plt.legend()
plt.grid(True)
plt.show()
