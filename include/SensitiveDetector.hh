#ifndef SENSITIVEDETECTOR_HH
#define SENSITIVEDETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4THitsCollection.hh"
#include "Hit.hh" // Include the Hit class

/**
 * @namespace G4Sim
 * @brief A namespace for fast simulation classes in Geant4.
/*/
namespace G4Sim {

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
    SensitiveDetector(const G4String& name, const G4String& volumeName);
    virtual ~SensitiveDetector();

    virtual void Initialize(G4HCofThisEvent* hce) override;
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    virtual void EndOfEvent(G4HCofThisEvent* hce) override;

    void AddCollection(G4String volumeName);

    G4double GetTotalEnergyDeposit() const { return fTotalEnergyDeposit; }
    G4String GetTag() const { return fTag; }

        // Add this getter function
    const std::vector<HitsCollection*>& GetHitsCollections() const {
        return fHitsCollections;
    }


private:
    //HitsCollection* fHitsCollection;
    std::vector<HitsCollection*> fHitsCollections;  // A vector to store multiple hits collections

    //G4THitsCollection<Hit>* fHitsCollection;
    G4double fTotalEnergyDeposit;
    G4int fHitsCollectionID;
    G4String fTag;
};

}

#endif
