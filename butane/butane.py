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

butane = read('butane.xyz')
butane.calc = NWChem(xc='PBE')
opt = BFGS(butane, trajectory="optimization.traj")
opt.run(fmax=0.01)
write('butane_opt.xyz', butane)
butane.get_potential_energy()