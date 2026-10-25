# PyNeut vs OpenMC validation suite

These scripts compare PyNeut against [OpenMC](https://docs.openmc.org).
Every OpenMC reference is generated **live**, from the same ENDF/B HDF5 files
PyNeut reads, with identical geometry and source. The numbers can therefore be
reproduced on any machine that has `openmc` and the data. The results and how
to read them are in the [validation report](../../docs/index.md).

**Requirements:**
- `pip install -e ".[validation]"`, which provides `openmc` and `matplotlib`.
- The ENDF/B-VIII.0 HDF5 files in `endfb/neutron/`.

The OpenMC library file `endfb/cross_sections.xml` is built automatically from
those files the first time a live run needs it.

## Scripts

| Script | What it checks |
|--------|----------------|
| `run_all.py` | **Main runner.** All 13 fixed-source cases plus the cross-section check. It prints the summary table and writes `validation_results.csv` and `validation_manifest.json`. `--plots` writes one spectrum panel per case to `docs/images/cases/`. |
| `validate_<iso>.py` | Case definitions for Pb208, Fe56, Be9, Al27, C12 and O16. Each file also runs on its own. |
| `validate_thermal.py` | The thermal and epithermal cases (B10, Cd112, C12). |
| `validate_xs.py` | The cross-section reader against `openmc.data.IncidentNeutron`, channel by channel. It needs no `cross_sections.xml`. |
| `validate_secondary.py` | (n,2n) emission isolated from transport: the energy and angle of the emitted neutron at a fixed incident energy, against OpenMC's parsed distribution. |
| `validate_mixture.py` | Per-collision isotope selection in a B10 + C12 sphere, with identical number densities in both codes. |
| `validate_keff.py` | k-eigenvalue: by default a sweep of bare U235 sphere radii. `--benchmarks` runs Godiva, Jezebel-23 and Jezebel against the published k_eff and needs U233/U234/U238 and Pu239–241/Ga data. |
| `validate_pincell.py` | k∞ for the BEAVRS hot-zero-power pin cell in an all-reflective cube. The matched-physics comparison uses free-gas H at 294 K and no S(α,β), over several seeds per code. With OpenMC's thermal-scattering data available, it also reports a realistic OpenMC value at 600 K with `c_H_in_H2O`. `--save-reference` rewrites `pincell_reference.json`. Needs the 36 pin-cell nuclides. |
| `run_all_vector.py` | The 13 fixed-source cases on the vector engine, graded against the committed OpenMC references (OpenMC is not re-run), with a scalar-vs-vector speed-up column. Writes `validation_results_vector.csv`. |
| `benchmark_pincell_vector.py` | The pin cell on the vector engine, compared against `pincell_reference.json`, plus throughput. |
| `benchmark_pincell_runtime.py` | Wall-clock time of PyNeut (scalar) vs OpenMC on the pin cell, with the same budget and core count. |
| `main_benchmark.py` | A minimal PyNeut-only Pb208 sphere timing run. |
| `_common.py` | Shared code: building the data library, running each code, batch statistics, the spectrum χ²/dof, grading and printing. |
| `../analytic_benchmarks.py` | PyNeut against exact results: Beer–Lambert attenuation and the elastic [α,1] / ξ checks. It needs data but not OpenMC. |

## Committed data

| File | Contents |
|------|----------|
| `validation_results.csv` | the reference fixed-source table, N = 10 000 |
| `validation_manifest.json` | a fingerprint of the ENDF snapshot `validation_results.csv` was produced with |
| `validation_results_vector.csv` | the vector-engine table, N = 100 000 |
| `pincell_reference.json` | 5-seed PyNeut and OpenMC k∞ for the matched-physics pin cell, with config and data note |
| `omc_pincell_seeds.json` | an earlier 6-seed OpenMC pin-cell reference (2026-07-07 data snapshot), kept for the paper |
| `openmc_runs/omc_*/` | `model.xml` and the final statepoint of each OpenMC run behind the tables. The scripts overwrite these on every run. |

Results are tied to one data snapshot. The same lattice can move by hundreds
of pcm across an evaluation refresh, so regenerate the references whenever the
library changes. A subset run (`--isotopes`) writes `*_subset` files, which
are git-ignored, and leaves the reference table alone.

## Running

```bash
cd Validation/OpenMC_Comparison

python run_all.py                                 # full suite, N = 10 000
python run_all.py --n 2000 --isotopes Pb208,Fe56  # quick subset
python run_all.py --plots                         # + per-case spectrum panels
python validate_xs.py                             # reader check only
python validate_secondary.py                      # (n,2n) emission diagnostic
python validate_mixture.py                        # B10+C12, default --n 40000
python validate_keff.py --n 2000 --gens 120 --inactive 20
python validate_keff.py --benchmarks              # Godiva, Jezebel-23, Jezebel
python validate_pincell.py --quick                # smoke test
python validate_pincell.py                        # 5 seeds/code, N = 3000 (slow)
python run_all_vector.py                          # vector engine, N = 100 000
python ../analytic_benchmarks.py
```

The pin cell is the slowest problem. Its reflective lattice has no leakage, so
thermal histories run for hundreds of collisions, and the scalar kernel is
about 100× heavier there than on the fast-metal spheres.
`benchmark_pincell_vector.py` runs the same problem on the vector engine.

## What is compared

1. **Integral metrics.** Leakage fraction and mean escape energy, each with a
   standard error: batch means for PyNeut, OpenMC's own batch statistics for
   OpenMC, with ratio error propagation for the energy. The energy uses the
   **same histogram** in both codes (`E_NBINS` bins from `energy_grid()`), so
   the binning bias cancels.
2. **Spectrum shape.** The whole leakage spectrum is compared with a
   coarse-group χ²/dof (`compare_spectra`): GOOD ≤ 2, FAIR ≤ 5, POOR > 5.
   This is reported separately from the integral status. PyNeut's spectrum
   is scored per leaking neutron (`pyneut.tallies`), as OpenMC's is.
   Results committed before that change binned each history's average
   escape energy instead, which distorts the shape whenever (n,xn) makes one
   history leak several neutrons.
3. **Fallback diagnostics.** The kernel counts which ENDF law sampled each
   (n,xn) secondary energy. The "Max%" column is the share that fell back to
   Maxwellian evaporation.

### Grading

```
z = |OpenMC - PyNeut| / sqrt(sem_OMC^2 + sem_PyN^2)
```

| Status | Condition (the worse of leakage and energy) |
|--------|---------------------------------------------|
| OK | z ≤ 2 |
| WARNING | 2 < z ≤ 4 |
| FAIL | z > 4 |
| NO_SIGMA | uncertainty unavailable (e.g. a single batch) |
