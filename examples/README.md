# Example inputs

Every case PyNeut is validated on, and one file per feature of the
[input format](../docs/input-format.md). Run any of them from the repository
root, with the ENDF/B-VIII.0 files in `endfb/neutron/` (or set
`PYNEUT_DATA`):

```bash
pyneut check examples/validation/eigenvalue/godiva.toml
pyneut run   examples/validation/eigenvalue/godiva.toml
```

The test suite loads, checks and runs every file here, against a synthetic
data library, so the examples stay in step with the code.

## Validation cases

These are the models behind the [validation report](../docs/index.md), with
the same geometry, densities, sources and spectrum bins as
`Validation/OpenMC_Comparison/`. The validation scripts also run OpenMC, so
use them for the comparison itself. These files let you run PyNeut's side of
any case on its own.

| File | Case |
|------|------|
| `validation/fixed_source/pb208_elastic.toml` | Pb208 sphere R=10, 2 MeV |
| `validation/fixed_source/pb208_5mev.toml` | Pb208 sphere R=10, 5 MeV |
| `validation/fixed_source/pb208_n2n.toml` | Pb208 sphere R=10, 14 MeV, (n,2n) |
| `validation/fixed_source/fe56_slab_14mev.toml` | Fe56 slab 4 cm thick, 14 MeV |
| `validation/fixed_source/fe56_sphere_1mev.toml` | Fe56 sphere R=10, 1 MeV |
| `validation/fixed_source/fe56_sphere_25kev.toml` | Fe56 sphere R=15, 25 keV |
| `validation/fixed_source/be9_n2n_14mev.toml` | Be9 sphere R=10, 14 MeV, (n,2n) |
| `validation/fixed_source/al27_slab_1mev.toml` | Al27 slab 1 cm thick, 1 MeV |
| `validation/fixed_source/c12_cyl_2mev.toml` | C12 cylinder R=10 × 20, 2 MeV |
| `validation/fixed_source/o16_sphere_10mev.toml` | O16 gas sphere R=20, 10 MeV |
| `validation/fixed_source/b10_epithermal_1kev.toml` | B10 sphere R=1, 1 keV |
| `validation/fixed_source/cd112_thermal.toml` | Cd112 sphere R=2, thermal |
| `validation/fixed_source/c12_thermal.toml` | C12 sphere R=5, thermal |
| `validation/mixture_b4c.toml` | B10 + C12 mixture sphere |
| `validation/eigenvalue/u235_sphere_r6.toml`, `_r8_7`, `_r11` | bare U235 metal spheres |
| `validation/eigenvalue/godiva.toml` | Godiva, HEU-MET-FAST-001 |
| `validation/eigenvalue/jezebel23.toml` | Jezebel-23, U233-MET-FAST-001 |
| `validation/eigenvalue/jezebel.toml` | Jezebel, PU-MET-FAST-001 |
| `validation/pincell_beavrs.toml` | BEAVRS pin cell, k∞ (slow) |

The validation reference table was produced with the scalar engine. To run a
case the same way, add `--engine scalar`.

## Features

| File | Shows |
|------|-------|
| `features/slab_detector.toml` | beam source, box surfaces, cells carved with `~` and `\|`, analog transport, detector tally |
| `features/mesh_b10_cylinder.toml` | capped cylinder, mesh absorption tally written as VTK (the `main.py` demo) |
| `features/infinite_medium.toml` | reflective boundaries |
| `features/region_operations.toml` | every region operator: `&`, `\|`, `~`, parentheses |
| `features/layered_shield.toml` | several materials in concentric shells, log-binned leakage spectrum |
| `features/priority_nesting.toml` | the same shield built from overlapping cells resolved by `priority` |
| `features/source_box_uniform.toml` | volume source in a box, uniform energy |
| `features/source_discrete_lines.toml` | discrete energy lines (DD + DT) |
| `features/source_maxwell_cell.toml` | source inside one cell, Maxwellian energy |
| `features/weight_fractions.toml` | material given by weight fractions |
| `features/analog_transport.toml` | analog instead of implicit-capture transport |
| `features/general_plane.toml` | tilted planes |
| `features/water_moderator.toml` | compound material by atom fractions, slowing down to thermal |
| `features/reflected_core.toml` | eigenvalue with a reflector, on the scalar engine, all cores |
