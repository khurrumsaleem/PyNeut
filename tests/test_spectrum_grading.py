"""The validation suite's spectrum grader (Validation/OpenMC_Comparison).

It needs neither nuclear data nor OpenMC: the spectra are made up.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "Validation", "OpenMC_Comparison"))

from _common import compare_spectra  # noqa: E402


# ------------------------------------------------------------------- grader
def _spectra(omc, pyn, rel=0.02, n=None):
    omc, pyn = np.asarray(omc, float), np.asarray(pyn, float)
    out = {"spectrum": omc, "spectrum_sem": rel * omc}
    if n:
        out["n_particles"] = n
    return out, {"spectrum": pyn, "spectrum_sem": rel * pyn}


def test_an_empty_analog_tail_is_not_a_disagreement():
    # PyNeut (implicit capture) resolves a tail worth ~0.9 analog counts in
    # 10 000 histories; OpenMC (analog) scored none of it, so its batch error
    # there is 0. That used to read as chi2/dof = 15 for B10 at 1 keV.
    body = np.linspace(1e-3, 5e-2, 12)
    tail = np.array([5e-9, 3e-8, 1e-7, 7e-7, 3e-6, 2e-5, 3e-5, 3e-5])
    omc, pyn = _spectra(np.r_[np.zeros(8), body], np.r_[tail, body])
    omc["n_particles"] = 10_000
    _, _, per_dof = compare_spectra(omc, pyn)
    assert per_dof < 1.0
    del omc["n_particles"]                  # the old scoring, for contrast
    assert compare_spectra(omc, pyn)[2] > 5


def test_a_real_shape_difference_is_still_caught():
    # 200-800 analog counts per group, so a +-20 % tilt is several sigma
    body = np.linspace(2e-2, 8e-2, 20)
    omc, pyn = _spectra(body, body * np.linspace(0.8, 1.2, 20), n=10_000)
    assert compare_spectra(omc, pyn)[2] > 5


def test_identical_well_populated_spectra_score_zero():
    body = np.linspace(1e-2, 5e-2, 20)
    omc, pyn = _spectra(body, body, n=10_000)
    chi2, ndof, _ = compare_spectra(omc, pyn)
    assert chi2 == 0.0 and ndof == 20
