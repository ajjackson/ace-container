from typing import Sequence

from julia.api import Julia

# This seems to be required or the Julia connection fails when HAL attempts it
jl = Julia(compiled_modules=False)  # flake8: E402 noqa

from ACEHAL.HAL import HAL
from ase import Atoms
from ase.calculators.emt import EMT
import ase.io
from sklearn.linear_model import BayesianRidge

def tweak_metadata(training_data: Sequence[Atoms]) -> None:
    for atoms in training_data:
        atoms.info["energy"] = atoms.get_potential_energy()
        atoms.arrays["forces"] = atoms.get_forces()
        atoms.info["virial"] = -atoms.get_stress(voigt=False) * atoms.cell.volume


def get_e0s(atoms, calc):
    elements = list(set(atoms.get_chemical_symbols()))

    e0s = {}

    for element in elements:
        at = Atoms(element, pbc=False)
        at.calc = calc
        e0s.update({element: at.get_potential_energy(force_consistent=True)})

        return e0s


def main():
    fit_configs = ase.io.read("rattled_ev.extxyz", index=":")
    supercells = ase.io.read("supercells.extxyz", index=":")
    tweak_metadata(fit_configs)

    data_keys = {"E": "energy", "F": "forces", "V": "virial", "Fmax": 15.0}

    E0s = get_e0s(fit_configs[0], EMT())
    weights = {"E_per_atom": 1e2, "F": 1e1, "V_per_atom": 1e1}

    solver = BayesianRidge(fit_intercept=True, compute_score=True)

    fixed_basis_info = {
        "elements": list(E0s),
        "cor_order": 2,
        "r_cut": 3.0,
        "smoothness_prior": ("algebraic", 2),
    }

    optimize_params = {"maxdeg": ("int", (3, 12))}

    basis_optim_kwargs = {
        "n_trials": 10,  # max number of basis opt iterations
        "timeout": 10000,
        "max_basis_len": 3000,
        "fixed_basis_info": fixed_basis_info,
        "optimize_params": optimize_params,
    }

    HAL(
        fit_configs,  # inital fitting database
        supercells,  # initial starting points)
        None,  # use ACE1x defaults
        solver,  # sklearn solver
        fit_kwargs={"E0s": E0s, "data_keys": data_keys, "weights": weights},
        n_iters=200,
        traj_len=10000,  # Max ML steps before evaluating calculator
        tol=0.2,  # relative uncertainty tolerance [0.2-0.4]
        tol_eps=0.2,  # regularising fraction of uncertainty [0.1-0.2]
        tau_rel=0.2,  # biasing strength [0.1-0.3]
        ref_calc=EMT(),
        dt_fs=2.0,  # MD timestep
        T_K=800,  # MD temperature
        T_timescale_fs=50,  # Langevin thermostate parameter
        P_GPa=0.0,  # Pressure
        swap_step_interval=0,  # atom swap MC step interval
        cell_step_interval=50,  # cell shape MC step interval
        basis_optim_kwargs=basis_optim_kwargs,
        basis_optim_interval=5,
        file_root="Al",
        test_fraction=0.1,
        traj_interval=5,
    )


if __name__ == "__main__":
    main()
