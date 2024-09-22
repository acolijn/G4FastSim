from setuptools import setup, find_packages

setup(
    name='G4Sim',
    version='0.1',
    packages=find_packages(where='analysis'),
    package_dir={'': 'analysis'},  # This tells it that the packages are in src/
)