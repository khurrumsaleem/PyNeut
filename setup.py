from setuptools import setup

setup(
    name="pyneut",
    version="0.1.0",
    description="Continuous-energy Monte Carlo neutron transport in Python, validated against OpenMC",
    url="https://github.com/ybadr16/PyNeut",
    license="MIT",
    packages=["pyneut"],
    python_requires=">=3.10",
    install_requires=["numpy", "h5py"],
    extras_require={
        "test": ["pytest", "trimesh", "rtree"],
        "validation": ["openmc", "matplotlib"],
    },
    entry_points={"console_scripts": ["pyneut = pyneut.cli:main"]},
)
