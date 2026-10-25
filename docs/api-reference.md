---
title: API Reference
layout: default
---

# API reference

The public interface of the `pyneut` package, which you import as
`from pyneut.<module> import ...`. Units are eV, cm, g/cm³ and atoms/cm³
throughout. The [Quick Start](QUICK_START.md) shows these pieces working
together.

**Modules:** [material](#material) · [medium](#medium) ·
[cross_section_read](#cross_section_read) · [settings](#settings) ·
[simulation](#simulation) · [tally](#tally) · [statistics](#statistics) ·
[mesh](#mesh) · [criticality](#criticality) ·
[vector_transport](#vector_transport) · [physics](#physics) ·
[vt_calc](#vt_calc) · [geometry](#geometry) ·
[random_number_generator](#random_number_generator)

---

## material

### `Material(name, density, atomic_mass, atomic_weight_ratio=None)`

A single-isotope material.

| Argument | Meaning |
|----------|---------|
| `density` | g/cm³ |
| `atomic_mass` | molar mass, g/mol (used for the number density and the free-gas mass) |
| `atomic_weight_ratio` | AWR = target mass / neutron mass (used for the collision kinematics) |

The constructor computes two attributes: `number_density` = ρ·N_A / M in
atoms/cm³, and `kg_mass`, the mass of one atom in kg. Pass `kg_mass` to
`VelocitySampler`.

```python
from pyneut.material import Material, NEUTRON_MASS_AMU
AWR = 206.190                                               # Pb208, from the data file
lead = Material("Pb208", 11.35, AWR * NEUTRON_MASS_AMU, atomic_weight_ratio=AWR)
```

### `Material.mixture(name, density, components, temperature=294)` → `list[Nuclide]`

This classmethod builds a multi-isotope composition from a mass density and
atom fractions. Each component is `(element, awr, atom_fraction)` or
`(element, awr, atom_fraction, molar_mass)`. The fractions don't need to be
normalised. When the molar mass is omitted it defaults to
`awr × NEUTRON_MASS_AMU`. Pass the result to `Region(composition=...)`.

### `Nuclide(element, number_density, atomic_weight_ratio, atomic_mass=None, temperature=294)`

One isotope inside a region's composition:

- `element` is the data key, such as `"U235"`, matching
  `endfb/neutron/U235.h5`.
- `number_density` is in atoms/cm³. Multiply atoms/b-cm by 10²⁴.
- The nuclide builds its own free-gas `sampler` at `temperature`.

### Constants

`NEUTRON_MASS_AMU = 1.00866491588` and `AVOGADRO = 6.022e23`.

---

## medium

Constructive solid geometry. A **surface** divides space into an inside,
where `evaluate(x, y, z) ≤ 0`, and an outside. A **Region** combines
surfaces, or other Regions, with a boolean operation.

### `Region(surfaces=None, operation="intersection", name=None, priority=0, is_void=False, element=None, composition=None)`

| Argument | Meaning |
|----------|---------|
| `surfaces` | list of surfaces and/or nested `Region`s |
| `operation` | `"intersection"` (all inside), `"union"` (any inside), `"complement"` (outside the first entry), or `"difference"` (inside the first, outside the second) |
| `priority` | where regions overlap, the highest priority wins |
| `is_void` | no interactions; neutrons stream to the next surface |
| `element` | single-isotope data key, which uses the `A`, `N` and `sampler` passed to the transport call |
| `composition` | `list[Nuclide]`. A mixture that overrides `element` |

A point inside no region is outside the geometry. A neutron that gets there
has leaked.

Methods: `contains(x, y, z, tolerance=1e-9)` → `bool`, and
`add_surface(surface)`.

### Surfaces

Every surface takes `boundary_type="transmission"` (the default: the neutron
crosses) or `"reflective"` (it is mirror-reflected). Each one has
`evaluate(x, y, z)`, `nearest_surface_method(x, y, z, u, v, w)` (the distance
along the ray, or `None`) and `normal(x, y, z)`.

| Class | Inside means |
|-------|--------------|
| `Plane(A, B, C, D)` | `A·x + B·y + C·z ≤ D`. The normal is normalised internally. `Plane(0, 0, 1, 5)` is z ≤ 5 and `Plane(-1, 0, 0, 0)` is x ≥ 0. |
| `Sphere(center, radius)` | distance from `center` ≤ `radius` |
| `Cylinder(axis, radius, center)` | infinite cylinder along `"x"`, `"y"` or `"z"` through `center`. Cap it with two planes. |
| `Box(x_min, x_max, y_min, y_max, z_min, z_max)` | the axis-aligned box. `Box` is a `Region` made of six planes, so you can put it in another region's `surfaces`. `boundary_type` applies to all six faces. |

```python
from pyneut.medium import Region, Cylinder, Plane, Box
shield = Region(surfaces=[Cylinder("z", 10, (0, 0, 0)),
                          Plane(0, 0, 1, 10), Plane(0, 0, -1, 10)],   # |z| <= 10
                name="Shield", priority=1, element="Pb208")
world = Region(surfaces=[Box(-20, 20, -20, 20, -20, 20)], name="Void", is_void=True)
```

---

## cross_section_read

### `CrossSectionReader(base_path, temperature="294K")`

This reads OpenMC-format HDF5 files from `<base_path>/neutron/<element>.h5`
and caches every table it loads. If the data has no `"294K"` block, it falls
back to `"293.6K"`.

> **Missing data is not an error.** An unknown element or MT, a missing file,
> or an energy below threshold all return 0. A misspelt nuclide therefore acts
> like vacuum. The aliases `C`/`Graphite` → C0, `Be` → Be9, `Al` → Al27,
> `Fe` → Fe56 and `Pb` → Pb208 each pick a single file. They don't expand an
> element into its natural isotopes.

| Method | Returns |
|--------|---------|
| `get_cross_section(element, mt, energy)` | microscopic σ in barns (lin–lin interpolation) |
| `get_cross_sections(element, energy, sampler, N, A, rng=None)` | macroscopic `(Σ_el, Σ_in, Σ_abs, Σ_fis, Σ_t)` in cm⁻¹. Σ_in sums MT 16, 17 and 51–91, and Σ_abs sums MT 102–107. With an `rng`, the unresolved-resonance probability tables are sampled. |
| `calculate_macroscopic_xs(mic, N)` | `mic · N · 1e-24`, or 0 for non-positive input |
| `get_inelastic_components(element, energy, N)` | `(Σ_in_total, [(mt, Σ, Q), …])` |
| `get_elastic_mu(element, energy, rng)` | a CM cosine from the tabulated MT 2 angular distribution, or an isotropic one if there is none |
| `get_secondary_correlated_sample(element, mt, energy, rng)` | `(E_out, mu, ok)` from Law 61 or Kalbach–Mann |
| `get_secondary_energy(element, mt, energy, rng)` | `E_out` from uncorrelated (Law 4) data, or `None` |
| `get_secondary_frame(element, mt)` | `"cm"` or `"lab"` |
| `get_nu(element, energy)` | total ν̄, or 0 if not fissile |
| `is_fissionable(element)` | `bool` |
| `get_cross_section_data(element, mt)` | the cached `{'energy', 'xs', 'threshold_energy', 'q_value'}`, or `None` |

---

## settings

### `Settings(mode="shielding", particles=1000, batches=10)`

| `mode` | `use_implicit_capture` | `weight_cutoff` | `roulette_survival_prob` | Use for |
|--------|:----------------------:|:---------------:|:------------------------:|---------|
| `"shielding"` | True | 1e-4 | 0.1 | leakage, absorption maps |
| `"criticality"` | False (analog) | 0 | 1.0 | k-eigenvalue, unweighted event counts |

The attributes are `mode`, `num_particles`, `batches` and the three columns
above. Any other mode raises `ValueError`.

---

## simulation

The per-history (scalar) engine. It is the reference implementation that the
vector engine is checked against.

### `simulate_single_particle(args)` → `dict`

Transports one source neutron and every secondary it produces through (n,xn).
In analog mode it also banks fission neutrons for the next generation. Its
single argument is a 10-tuple, so it can be passed straight to `Pool.map`:

| # | Item | Notes |
|---|------|-------|
| 1 | `state` | `dict`: `x, y, z` (cm); `theta` (polar angle from +z) and `phi` (azimuth from +x), in rad; `energy` (eV); `weight` (normally 1.0); `has_interacted` (`False`) |
| 2 | `reader` | `CrossSectionReader` |
| 3 | `mediums` | `list[Region]` |
| 4 | `A` | AWR of the `element=` regions |
| 5 | `N` | number density of the `element=` regions (atoms/cm³) |
| 6 | `sampler` | `VelocitySampler` of the `element=` regions |
| 7 | `region_bounds` | detector box `(x_min, x_max, y_min, y_max, z_min, z_max)`, or `None`. It is tested only where the neutron stops, at collisions and crossings. |
| 8 | `track_coordinates` | `bool`: record the trajectory |
| 9 | `rng` | `RNGHandler` for this history |
| 10 | `settings` | `Settings` |

Items 4–6 are ignored for regions that have a `composition`.

| Result key | Meaning |
|------------|---------|
| `absorbed_weight` | total weight removed, by implicit capture or analog absorption |
| `absorbed_coords` | `[(x, y, z, weight_removed), …]`, for mesh scoring |
| `final_weight` | total weight that leaked |
| `final_energy` | leaked-weight-averaged escape energy (eV), or 0 if nothing leaked |
| `region_detected` | `bool`: the detector box was reached |
| `fissioned` | number of fission terminations |
| `fission_sites` | the next-generation source states (analog mode) |
| `physics_counts` | event counters, e.g. collisions, (n,xn) sampling law, `fission_production` |
| `source_uncollided` | the source neutron leaked without colliding |
| `trajectory` | `[(x, y, z), …]` for the history and its descendants, or `None` |
| `result` | always `"simulated"` |

`run_particle_kernel(state, reader, mediums, A, N, sampler, region_bounds, track_coordinates, rng, settings, counters=None, fission_bank=None)`
transports one neutron (no descendants) and returns
`(status, children, trajectory, absorbed_weight, absorbed_coords)`.

`THERMAL_KINEMATICS_CUTOFF = 10.0` eV. Below it, elastic scattering uses
free-gas kinematics.

---

## tally

### `Tally()`

This accumulates per-history results in source order, which the batch-means
errors depend on.

| Member | Meaning |
|--------|---------|
| `merge_partial_results(res)` | add one `simulate_single_particle` result |
| `statistics(num_particles, n_batches=10)` | `{'leakage', 'leakage_sem', 'avg_energy', 'avg_energy_sem'}` from batch means |
| `print_summary(num_particles, n_batches=10)` | prints the above, plus the absorbed weight and the detector count |
| `get_results()` | all the raw lists in one dict |
| `leaked_weights`, `energy_spectrum` | per-history `final_weight` and `final_energy` |
| `absorbed_coordinates` | every `absorbed_coords` entry |
| `region_count` | histories that reached the detector box |
| `results["absorbed"]` | total absorbed weight |

`results["escaped"]` and `results["killed"]` are counters for the legacy
`update()` path and stay at 0 under `merge_partial_results`. Use
`statistics()` for leakage.

---

## statistics

| Function | Purpose |
|----------|---------|
| `escape_statistics(leaked_weight, escape_energy, n_source, n_batches=10, e_bins=None)` | batch-means leakage and mean escape energy, each ± its standard error, from per-history arrays. With `e_bins`, the energy comes from a weighted histogram on the bin midpoints, the same estimator the OpenMC comparison uses. |
| `mean_sem(batch_values)` | `(mean, standard error of the mean)`, ignoring NaN batches |
| `batch_edges(n_items, n_batches)` | contiguous batch boundaries |
| `ratio_mean_uncertainty(values, weights, weight_sigmas)` | weighted mean, with the error propagated from independent weight errors |

---

## mesh

### `MeshTally(x_bounds, y_bounds, z_bounds, dims)`

A regular Cartesian mesh. `score(x, y, z, weight)` adds to the voxel that
contains the point, and points outside the mesh are dropped.
`write_vtk(filename)` writes legacy-VTK `STRUCTURED_POINTS` cell data, which
ParaView opens directly.

---

## criticality

### `run_keff(reader, mediums, initial_source, A, N, sampler, settings, generations=50, inactive=10, seed=1, pool=None, verbose=False)` → `dict`

Fission-source power iteration on the scalar engine. `settings.mode` must be
`"criticality"`. Each generation transports its source (in parallel when
`pool` is given), then resamples the fission bank back to the starting size.
k is averaged over the active generations. Per-generation progress is logged
to `pyneut.criticality`, at INFO with `verbose=True` and at DEBUG otherwise. To
see it, configure logging, for example with
`logging.basicConfig(level=logging.INFO)`. If the fission chain dies out, a
warning is logged and the iteration stops early.

| Key | Meaning |
|-----|---------|
| `k_eff`, `k_sem` | collision estimator w·Σ νΣ_f/Σ_t (primary) |
| `k_source`, `k_source_sem` | generation estimator: fission neutrons produced per neutron started |
| `active_k` | per-generation collision k over the active generations |
| `cycles` | `[(generation, k_collision, k_source, is_active), …]` |

### `uniform_sphere_source(center, radius, n, rng, element="U235", energy=None)`

Returns `n` isotropic source states spread uniformly through a sphere. Their
energies follow the element's Watt spectrum unless `energy` is given.

---

## vector_transport

The event-based engine. It moves the whole bank one event (a flight, then a
crossing or collision) per NumPy pass, with the same physics as `simulation`.

A **source** is a dict of equal-length arrays `x, y, z, u, v, w, energy`
(unit direction cosines), with an optional `weight`. `rng` is a
`numpy.random.Generator`. Material regions need a `composition`, or you can
pass `legacy=(A, N, sampler)` to back `element=` regions.

| Function | Returns |
|----------|---------|
| `run_transport(reader, mediums, source, rng, settings, legacy=None, fission_bank=None, max_events=1_000_000)` | `{'leaked_weight', 'escape_energy'}` (per-history arrays for `escape_statistics`), plus `'leakage'`, `'absorbed_weight'`, `'counters'` and `'events'`. In analog mode with a `fission_bank` list, fission neutrons are banked as array tuples. |
| `run_transport_parallel(base_path, mediums, source, settings, n_workers=None, seed=0, temperature="294K", legacy=None)` | the same dict. The bank is split into one contiguous chunk per worker, and each worker builds its own reader from `base_path`. |
| `run_keff_vector(reader, mediums, initial_source, settings, generations=50, inactive=10, seed=1)` | the same dict as `criticality.run_keff` |
| `run_keff_parallel(base_path, mediums, sources, settings, seeds, generations=50, inactive=10, temperature="294K", n_workers=None)` | `{'k_mean', 'k_spread_sem', 'per_seed_k', 'per_seed'}`: independent replicas, with the error bar taken from their spread |
| `run_streaming(reader, mediums, source, rng, legacy=None)` | uncollided transmission: any collision ends the history (a Beer–Lambert check) |

Building blocks (array versions of the scalar routines): `VectorXS`,
`nearest_crossings`, `resolve_media`, `contains_many`, `surface_distances`,
`surface_normals`, `rotate_directions`, `static_elastic`, `free_gas_elastic`,
`discrete_inelastic`, `sample_velocity_many`, `sample_mu_many`,
`sample_watt_many`.

---

## physics

Scalar kinematics, used by the per-history kernel and the tests.

| Function | Returns |
|----------|---------|
| `elastic_scattering(E, A, sampler, rng)` | `(E′, μ_cm, μ_lab)` from free-gas elastic scattering, isotropic in the CM frame |
| `inelastic_scattering(E, A, Q, sampler, rng)` | `(E′, μ_cm, μ_lab)` for a discrete level. It falls back to elastic below threshold. |
| `sample_new_direction_cosines(u, v, w, μ_lab, rng)` | `(u′, v′, w′, φ)` |
| `sample_maxwellian(T_eV, rng)` | an energy from √E·e^(−E/T) |
| `get_nuclear_temperature(E, A)` | the Fermi-gas T = √(U/a), with a = A/8 MeV⁻¹, in eV |
| `sample_watt_spectrum(rng, a, b)` | a Watt-spectrum energy (eV), with `a` in eV and `b` in 1/eV |
| `watt_params_for(element)` | `(a, b)` for U235, U233 or Pu239. Anything else defaults to U235. |
| `calculate_E_cm_prime`, `calculate_E_prime`, `get_vectors_and_vcm`, `scatter_isotropic_cm` | the lower-level pieces of the above |

---

## vt_calc

### `VelocitySampler(mass, temperature=294)`

The free-gas target sampler. `mass` is the mass of one target atom in kg
(`Material.kg_mass`). `sample_velocity(vn, rng, max_attempts=1000, return_mu=False)`
samples the target speed (m/s) for a neutron moving at speed `vn`. The `rng`
argument, the history's `RNGHandler`, is required, so every draw comes from
that history's reproducible stream. With
`return_mu=True` it returns `(v_t, μ)`, the speed and its cosine to the
neutron direction, which must be used together.

---

## geometry

| Function | Returns |
|----------|---------|
| `calculate_nearest_crossing(state, regions, u, v, w)` | `(point, region, distance, surface)` for the next crossing, or `(None, None, inf, None)` |
| `calculate_nearest_boundary(state, regions, u, v, w)` | the same, without `surface` |
| `reflect_direction(u, v, w, nx, ny, nz)` | the mirror-reflected direction |
| `get_primitive_surfaces(surfaces)` | the primitive surfaces of a nested surface list, flattened |
| `count_coordinates_in_boundary(coords, x_bounds, y_bounds, z_bounds)` | the number of points inside the box |

---

## random_number_generator

### `RNGHandler(seed=None)`

A thin wrapper around `numpy.random.default_rng`. Methods: `random()`,
`random_array(size)`, `uniform(low, high)`, `log_uniform(scale)` (an
exponential sample, −ln(1−ξ)/scale) and `choice(...)`. Give every history its
own `RNGHandler`, for example with seed `base + i`, so that results are
reproducible no matter how `Pool` schedules the work.
