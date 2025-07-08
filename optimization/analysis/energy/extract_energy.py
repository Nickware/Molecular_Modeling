import re
import matplotlib.pyplot as plt

energies_ha = [] # Energies in Hartree

with open("optimizacion.out", "r") as f:
    for line in f:
        if "Total DFT energy" in line:
            # Use regex to extract the first valid floating-point number
            match = re.search(r"[-+]?\d*\.\d+|\d+", line)
            if match:
                energy = float(match.group())
                energies_ha.append(energy)
            else:
                print(f"Warning: Could not extract energy from line: {line.strip()}")

# Convert Hartree to eV (optional, ASE uses eV by default)
energies_ev = [e * 27.2114 for e in energies_ha]

# Convert Hartree to eV (optional, ASE uses eV by default)
energies_ev = [e * 27.2114 for e in energies_ha]
print(energies_ev)
# Create plot
plt.figure(figsize=(8, 5))
plt.plot(energies_ev, 'bo-', markersize=5, label="DFT Energy (eV)")
# Customize plot
plt.xlabel("Optimization Step", fontsize=12)
plt.ylabel("Energy (eV)", fontsize=12)
plt.title("DFT Energy Convergence", fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()

plt.show()