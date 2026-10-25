# Contributing to PyNeut

Thanks for your interest in PyNeut. Contributions, bug reports and questions are
all welcome. PyNeut is a small, transparent Monte Carlo neutron transport code
intended for education and research prototyping, so readability and a clear
validation trail matter as much as new features.

## Reporting bugs and asking questions

- Please open an issue on the [GitHub issue tracker](https://github.com/ybadr16/PyNeut/issues).
- For a **bug**, include: what you ran (geometry/source/settings or a minimal
  script), what you expected, what happened, and your Python, NumPy and `h5py`
  versions. A minimal reproducible example is the fastest path to a fix.
- For a **question** or usage help, open an issue with the `question` label;
  there is no separate mailing list.

## Development setup

PyNeut needs **Python 3.10+**. CI tests on 3.10 and 3.12. At runtime it
depends only on NumPy and `h5py`.

```bash
git clone https://github.com/ybadr16/PyNeut.git
cd PyNeut
pip install -e ".[test]"          # numpy, h5py + pytest, trimesh, rtree
pip install -e ".[validation]"    # + openmc, matplotlib for Validation/OpenMC_Comparison/
```

The ENDF/B-VIII.0 nuclear data is **not** in the repository, because it is
large and publicly available. Download the OpenMC-format HDF5 files into
`endfb/neutron/` as the README describes before you run transport or
validation.

## Running the tests

From the repository root:

```bash
pytest -q
```

Without nuclear data, the 42 tests that read ENDF files skip and the rest run.
This is also what CI does. If you have the data, please run the full suite
before opening a pull request, since the skipped tests cover the transport
end to end. If you add
or change physics, add a test that pins the behaviour (a unit test for kinematics
or geometry, or an entry in the validation suite for a transport-level check).

If you add a feature to the input format, also add an input file under
`examples/` that uses it. `tests/test_input_format.py` loads, checks and runs
every example against a synthetic data library (`tests/synthetic_data.py`),
so CI covers the examples too. That library has flat, made-up cross sections,
so these tests check the code paths, not physics values.

## Submitting changes

1. Fork the repository and create a topic branch off `main`.
2. Keep changes focused and match the surrounding style. The scalar kernel
   in `simulation.py` favours small, readable functions over
   micro-optimisation.
3. Add or update tests, and run `pytest -q`.
4. If your change affects results validated against OpenMC, analytic benchmarks or
   the ICSBEP cases, please note the before/after numbers in the pull request.
5. Open a pull request describing **what** changed and **why**, and link any
   related issue.

## Scope

PyNeut is deliberately not a production code. The contributions that fit best
improve transparency, validation coverage, documentation or physics fidelity.
Performance work goes in the vector engine (`vector_transport.py`). The scalar
kernel stays the readable reference, and any faster path has to be checked
against it statistically, the way `tests/test_vector_transport.py` does. If in
doubt, open an issue to discuss it before investing significant effort.

## Code of conduct

Please be respectful and constructive in all project spaces. Harassment or
demeaning behaviour will not be tolerated.

## License

By contributing, you agree that your contributions will be licensed under the
[MIT License](LICENSE) that covers the project.
