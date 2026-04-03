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

propane = read('propane.xyz')
propane.calc = NWChem(xc='PBE')
opt = BFGS(propane, trajectory="optimization.traj")
opt.run(fmax=0.01)
write('propane_opt.xyz', propane)
propane.get_potential_energy()