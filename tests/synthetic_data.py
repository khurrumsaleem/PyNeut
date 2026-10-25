"""A tiny synthetic nuclear-data library for tests that run without ENDF data.

Each file has the OpenMC HDF5 layout the reader expects (energy grid, MT 2 /
102 / 18 cross sections, total nu for fissile nuclides) with flat, made-up
cross sections, and the nuclide's real atomic-weight ratio. Results computed
on it exercise the code paths, not the physics.
"""
import os

import h5py
import numpy as np

NEUTRON_MASS_AMU = 1.00866491588

#: atomic masses (u) of every nuclide the examples and tests use
MASSES = {
    "H1": 1.00782503, "H2": 2.01410178, "He3": 3.01602932, "He4": 4.00260325,
    "Be9": 9.0121831, "B10": 10.01293695, "B11": 11.00930536, "C12": 12.0,
    "O16": 15.99491462, "O17": 16.99913176, "O18": 17.99915961,
    "Al27": 26.98153841,
    "Cr50": 49.94604183, "Cr52": 51.94050623, "Cr53": 52.94064815, "Cr54": 53.93887916,
    "Fe54": 53.93960899, "Fe56": 55.93493633, "Fe57": 56.93539284, "Fe58": 57.93327443,
    "Mn55": 54.9380439, "Ni58": 57.9353429, "Ni60": 59.9307864,
    "Ga69": 68.9255735, "Ga71": 70.9247026,
    "Zr90": 89.9046977, "Zr91": 90.9056396, "Zr92": 91.9050347, "Zr94": 93.9063108,
    "Zr96": 95.9082714,
    "Cd112": 111.9027629,
    "Sn112": 111.9048239, "Sn114": 113.9027827, "Sn115": 114.9033447,
    "Sn116": 115.9017428, "Sn117": 116.902954, "Sn118": 117.9016066,
    "Sn119": 118.9033112, "Sn120": 119.9022016, "Sn122": 121.9034438,
    "Sn124": 123.9052766,
    "Pb207": 206.9758973, "Pb208": 207.9766525,
    "U233": 233.0396355, "U234": 234.0409523, "U235": 235.0439301,
    "U236": 236.0455682, "U238": 238.0507884,
    "Pu239": 239.0521634, "Pu240": 240.0538135, "Pu241": 241.0568515,
}

FISSILE = {"U233", "U235", "Pu239", "Pu241"}
ABSORBERS = {"B10": 3.0, "Cd112": 2.0}


def write_library(directory, nuclides=None):
    """Write <directory>/neutron/<nuclide>.h5 for `nuclides` (default: all)."""
    out = os.path.join(directory, "neutron")
    os.makedirs(out, exist_ok=True)
    energy = np.logspace(-5, np.log10(2e7), 200)
    for name in nuclides or MASSES:
        fissile = name in FISSILE
        xs = {2: 4.0, 102: ABSORBERS.get(name, 0.1 if fissile else 0.02)}
        if fissile:
            xs[18] = 1.3
        with h5py.File(os.path.join(out, f"{name}.h5"), "w") as f:
            f.attrs["filetype"] = np.bytes_("data_neutron")
            root = f.create_group(name)
            root.attrs["atomic_weight_ratio"] = MASSES[name] / NEUTRON_MASS_AMU
            root.create_dataset("energy/294K", data=energy)
            for mt, value in xs.items():
                d = root.create_dataset(f"reactions/reaction_{mt:03}/294K/xs",
                                        data=np.full(energy.size, value))
                d.attrs["threshold_idx"] = 0
            if fissile:
                root.create_dataset("total_nu/yield",
                                    data=np.vstack([energy, np.full(energy.size, 2.5)]))
    return directory
