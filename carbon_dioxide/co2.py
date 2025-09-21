# Structure optimization of CO2 molecule with NWChem calculator
from ase import Atoms
from ase.io import read
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write

co2 = read('co2.xyz')
#h2 = Atoms('co2',
#           positions=[[0, 0, 0],
#                      [0, 0, 0.7]])
co2.calc = NWChem(xc='PBE')
opt = BFGS(co2, trajectory="optimization.traj")
opt.run(fmax=0.01)
write('co2_opt.xyz', co2)
co2.get_potential_energy()