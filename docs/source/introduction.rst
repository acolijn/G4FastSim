Installation
============

To install the package, run the following command::
    
    git clone https://github.com/acolijn/G4Sim.git

Your environment should include Geant4 and ROOT installations.

To initialize the environment on teh Nikhef computing cluster::
    
    conda activate g4

Usage
=====

Run the code with the following command::

    cd <PROJECT_BASE_DIR>/run
    python run_simulation.py -json <JSON_FILE> -n <NUMBER_OF_EVENTS> <OPTIONAL_COMMANDS >

