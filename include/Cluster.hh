#ifndef CLUSTER_HH
#define CLUSTER_HH

#include <vector>
#include "G4ThreeVector.hh"
#include "Hit.hh"

/**
 * @namespace G4Sim
 * @brief A namespace for fast simulation classes in Geant4.
/*/
namespace G4Sim {

/**
 * @struct Cluster
 * @brief A struct representing a cluster of hits in the simulation.
 *
 * The Cluster struct represents a cluster of hits in the simulation, with a position, energy deposit, time, and a vector of hits.
 */
struct Cluster {
    G4ThreeVector position;
    G4double energyDeposit;
    G4double time;
    std::vector<Hit*> hits;
};

} // namespace G4Sim

#endif
