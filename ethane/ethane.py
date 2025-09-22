# Structure optimization of ethane molecule with NWChem calculator
# (This is a simplified version without detailed input parameters)
# Note: Ensure NWChem is properly installed and configured in your environment.
# You may need to adjust the calculator settings based on your NWChem installation.
# Import necessary modules

from ase import Atoms
from ase.io import read
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write

ethane = read('ethane.xyz')
ethane.calc = NWChem(xc='PBE')
opt = BFGS(ethane, trajectory="optimization.traj")
opt.run(fmax=0.01)
write('ethane_opt.xyz', ethane)
ethane.get_potential_energy()