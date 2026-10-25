from setuptools import setup

setup(
    name="pyneut",
    version="0.1",
    description="Continuous-energy Monte Carlo neutron transport in Python, validated against OpenMC",
    url="https://github.com/ybadr16/PyNeut",
    license="MIT",
    # The package is imported as `src` (e.g. `from src.simulation import ...`).
    packages=["src"],
    python_requires=">=3.10",
    install_requires=["numpy", "h5py"],
    extras_require={
        "test": ["pytest", "trimesh", "rtree"],
        "validation": ["openmc", "matplotlib"],
    },
)
