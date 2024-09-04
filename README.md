# G4FastSim
Variance reduction for gamma ray MC in Geant4

## Compilation:

mkdir build

cd build/

cmake -DCMAKE_PREFIX_PATH=<Geant4_installation_dir> ../

## Running:

cd run

python run_simulation.py <arguments>