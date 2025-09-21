# Structure optimization of methane molecule with NWChem calculator
from ase import Atoms
from ase.io import read
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write

methane = read('methane.xyz')
methane.calc = NWChem(xc='PBE')
opt = BFGS(methane, trajectory="optimization.traj")
opt.run(fmax=0.01)
write('methane_opt.xyz', methane)
methane.get_potential_energy()