# Structure optimization of CO2 molecule with Quantum ESPRESSO calculator
from ase import Atoms
from ase.io import read
from ase.optimize import BFGS
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.io import write

# Leer la estructura inicial
co2 = read('co2.xyz')

# Configurar los parámetros de Quantum ESPRESSO
pseudopotentials = {
    'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
    'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF'
}

# Optionally create profile to override paths in ASE configuration:
profile = EspressoProfile(
    command='/home/jntorresr/workspace/q-e/build/bin/pw.x', 
    pseudo_dir='/home/jntorresr/workspace/q-e/pseudo'
)

input_data = {
    'control': {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'prefix': 'co2',
        'outdir': './',
        'pseudo_dir': './pseudos/',  # Directorio donde están los pseudopotenciales
        'tprnfor': True,
        'tstress': True
    },
    'system': {
        'ecutwfc': 40,  # Energía de corte en Ry
        'ecutrho': 320,  # Energía de corte para densidad
        'input_dft': 'PBE',
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.01
    },
    'electrons': {
        'conv_thr': 1.0e-8,
        'mixing_beta': 0.7
    }
}
 
# Configurar el calculador de Quantum ESPRESSO
co2.calc = Espresso(profile=profile, pseudopotentials=pseudopotentials,kpts=(1, 1, 1),  # Puntos k para molécula aislada
    calculation='scf'
)
"""
# Optimización de estructura
opt = BFGS(co2, trajectory="optimization.traj")
opt.run(fmax=0.01)

# Guardar estructura optimizada
write('co2_opt.xyz', co2)

# Obtener energía potencial
energy = co2.get_potential_energy()
print(f"Energía potencial optimizada: {energy} eV") """