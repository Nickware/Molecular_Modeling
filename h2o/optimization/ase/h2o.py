# Structure optimization of water molecule with NWChem and ASE
# This script uses ASE to optimize the geometry of a water molecule (H2O)
# using the NWChem calculator. It sets up the water molecule with
# specified bond length and angle, performs the optimization,
# and writes the optimized structure to an XYZ file.
# The script also calculates the potential energy of the optimized structure.
# The script uses the BFGS algorithm for optimization and sets a force
# convergence criterion of 0.02 eV/Å.
# The script requires the ASE package and NWChem to be installed and properly configured.
# The script assumes that NWChem is available in the system path.
# The script uses the PBE functional for the calculation.
# The script uses the ASE Atoms class to create the water molecule
# and the ASE BFGS class for optimization

# The script uses the ASE NWChem calculator to perform the calculation. 
# The script uses the ASE numpy module to perform numerical calculations.
# The script uses the ASE io module to write the optimized structure to an XYZ file.
# The script uses the ASE optimize module to perform the optimization.
# The script uses the ASE calculators module to perform the calculation.
# The script uses the ASE atoms module to create the water molecule.

from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write
import numpy as np

d = 0.9575
t = np.pi / 180 * 104.51

water = Atoms('H2O',
              positions=[(d, 0, 0),
                         (d * np.cos(t), d * np.sin(t), 0),
                         (0, 0, 0)])

water.calc = NWChem(xc='PBE')
opt = BFGS(water, trajectory="optimization.traj")
opt.run(fmax=0.01)
write('H20.xyz', water)
water.get_potential_energy()
# The script uses the ASE io module to write the optimized structure to an XYZ file.
