import ase.build
from ase.calculators.emt import EMT
from ase.filters import UnitCellFilter
import ase.io
from ase.optimize.bfgs import BFGS
import numpy as np
                                    

ref_atoms = ase.build.bulk('Al', cubic=True)
ref_atoms.calc = EMT()
uc = UnitCellFilter(ref_atoms)
opt = BFGS(uc)
opt.run(fmax=1e-5)
ref_atoms.write('opt_cubic.extxyz')

for i, scale in enumerate(np.linspace(0.97, 1.03, 15)):
    atoms = ref_atoms.copy()
    atoms.set_cell(atoms.cell * scale, scale_atoms=True)
    atoms.rattle(seed=i)

    atoms.calc = EMT()
    forces = atoms.get_forces()
    virial = -atoms.get_stress(voigt=False) * atoms.cell.volume

    atoms.info["energy"] = atoms.get_potential_energy(force_consistent=True)
    atoms.info["virial"] = virial

    ase.io.write('rattled_ev.extxyz', atoms, append=True)

for i, supercell in enumerate([
                [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
                [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                [[2, 0, 0], [2, 2, 0], [0, 0, 2]],
                [[3, 0, 0], [2, 2, 0], [0, 0, 3]]
]):

        sc = ase.build.make_supercell(ref_atoms, supercell)
        sc.rattle(seed=i)
        sc.calc = EMT()

        forces = sc.get_forces()
        virial = -sc.get_stress(voigt=False) * sc.cell.volume

        sc.info["energy"] = sc.get_potential_energy(force_consistent=True)
        sc.info["virial"] = virial

        ase.io.write('supercells.extxyz', sc, append=True)
