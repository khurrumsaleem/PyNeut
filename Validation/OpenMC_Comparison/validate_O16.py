"""
PyNeut vs OpenMC — O16 validation cases

Cases
-----
  o16_sphere_10mev : O16 sphere R=20 cm, density=0.0014 g/cm3, E=10 MeV
  o16_sphere_14mev : O16 sphere R=15 cm, density=1.14 g/cm3 (liquid oxygen),
                     E=14 MeV (MT 22, 28 and the (n,a)/(n,p) partials)

OpenMC reference is computed live against the local data library.
"""
from _common import validate_case, print_result

_CASES = [
    {'name': 'o16_sphere_10mev', 'el': 'O16', 'rho': 0.0014, 'A': 15.85751,
     'E': 10e6, 'geo': 'sphere', 'dims': [20.0]},
    {'name': 'o16_sphere_14mev', 'el': 'O16', 'rho': 1.14, 'A': 15.85751,
     'E': 14e6, 'geo': 'sphere', 'dims': [15.0]},
]


def run(n_particles=10_000):
    results = []
    for case in _CASES:
        print(f"  Running {case['name']} …")
        r = validate_case(case, n_particles, isotope='O16')
        print_result(r)
        results.append(r)
    return results


if __name__ == '__main__':
    print('\n=== O16 Validation ===')
    run()
