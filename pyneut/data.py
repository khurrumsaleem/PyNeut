"""Locating nuclear data and reading per-nuclide constants from it."""
import glob
import hashlib
import os
import re

import h5py

#: Environment variable naming the data directory when an input gives none.
DATA_ENV = "PYNEUT_DATA"

_NUCLIDE_RE = re.compile(r"^[A-Z][a-z]?\d{1,3}(_m\d)?$")


def resolve_data_dir(explicit=None, base_dir="."):
    """The data directory for a run.

    An explicit path is taken relative to `base_dir` (the input file's
    folder). Otherwise $PYNEUT_DATA, otherwise ./endfb in the working
    directory.
    """
    if explicit:
        path = explicit if os.path.isabs(explicit) else os.path.join(base_dir, explicit)
    else:
        path = os.environ.get(DATA_ENV) or "endfb"
    return os.path.abspath(path)


def nuclide_file(data_dir, name):
    return os.path.join(data_dir, "neutron", f"{name}.h5")


def check_nuclide_name(name):
    """Return None if `name` looks like a nuclide (U235, Am242_m1), else a message."""
    if _NUCLIDE_RE.match(name):
        return None
    if re.fullmatch(r"[A-Z][a-z]?", name):
        return (f"{name!r} is an element; natural elements are not expanded, "
                "so list its isotopes by mass number (e.g. Fe56)")
    return f"{name!r} is not a nuclide name like 'U235' or 'Pb208'"


def read_awr(data_dir, name):
    """Atomic-weight ratio of `name` from its OpenMC HDF5 file.

    Raises FileNotFoundError if the file is missing and KeyError if it has no
    `atomic_weight_ratio` attribute.
    """
    path = nuclide_file(data_dir, name)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    with h5py.File(path, "r") as f:
        if f"neutron/{name}" in f:
            group = f[f"neutron/{name}"]
        elif name in f:
            group = f[name]
        else:
            roots = [k for k in f.keys() if k != "filetype"]
            group = f[roots[0]] if roots else f
        if "atomic_weight_ratio" not in group.attrs:
            raise KeyError(f"{path} has no 'atomic_weight_ratio' attribute")
        return float(group.attrs["atomic_weight_ratio"])


def data_manifest(data_dir):
    """Fingerprint of a nuclear-data directory.

    A result is a statement about a code AND a snapshot of the evaluation
    files, so every recorded result carries this digest. It is a SHA-256 over
    the sorted (name, size) pairs of <data_dir>/neutron/*.h5: cheap, and
    sensitive to any file being added, removed or replaced. Returns
    {'n_files', 'total_bytes', 'sha256'} (digest truncated to 16 hex digits).
    """
    files = sorted(glob.glob(os.path.join(data_dir, "neutron", "*.h5")))
    h = hashlib.sha256()
    total = 0
    for f in files:
        size = os.path.getsize(f)
        total += size
        h.update(f"{os.path.basename(f)}:{size}\n".encode())
    return {"n_files": len(files), "total_bytes": total,
            "sha256": h.hexdigest()[:16]}
