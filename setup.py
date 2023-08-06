from setuptools import setup, find_packages

setup(
    name="moleculegraph",
    version="0.0.1",
    description="A graph representation for molecules and everything else for data mapping",
    url="https://github......",
    author="Maximilian Fleck",
    author_email="fleck@itt.uni-stuttgart.de",
    license="BSD 2-clause",
    packages=find_packages(),
    install_requires=["ase", "numpy", "networkx", "rapidfuzz", "toml", "rdkit"],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.8",
    ],
)
