#ifndef SENSITIVEDETECTOR_HH
#define SENSITIVEDETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4THitsCollection.hh"
#include "Hit.hh" // Include the Hit class

/**
 * @namespace G4FastSim
 * @brief A namespace for fast simulation classes in Geant4.
/*/
namespace G4FastSim {

/**
 * @class SensitiveDetector
 * 
 * @brief This class represents a sensitive detector in the Geant4 simulation.
 * 
 * The SensitiveDetector class is derived from the G4VSensitiveDetector class and is responsible for
 * detecting and processing hits in the simulation. It provides methods for initializing the hits collection,
 * processing hits, and ending the event. It also includes a method to retrieve the total energy deposit.
 * 
 * @note This class assumes the existence of a HitsCollection class and a Hit class.
 */
class SensitiveDetector : public G4VSensitiveDetector {
public:
    SensitiveDetector(const G4String& name, const G4String& hitsCollectionName);
    virtual ~SensitiveDetector();

    virtual void Initialize(G4HCofThisEvent* hce) override;
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    virtual void EndOfEvent(G4HCofThisEvent* hce) override;

    G4double GetTotalEnergyDeposit() const { return fTotalEnergyDeposit; }

private:
    HitsCollection* fHitsCollection;
    //G4THitsCollection<Hit>* fHitsCollection;
    G4double fTotalEnergyDeposit;
    G4int fHitsCollectionID;
};

}

#endif
