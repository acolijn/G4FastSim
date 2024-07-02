#ifndef SENSITIVEDETECTOR_HH
#define SENSITIVEDETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4THitsCollection.hh"
#include "Hit.hh" // Include the Hit class

namespace G4FastSim {

class SensitiveDetector : public G4VSensitiveDetector {
public:
    SensitiveDetector(const G4String& name);
    virtual ~SensitiveDetector();

    virtual void Initialize(G4HCofThisEvent* hce) override;
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    virtual void EndOfEvent(G4HCofThisEvent* hce) override;

    G4double GetTotalEnergyDeposit() const { return fTotalEnergyDeposit; }

private:
    HitsCollection* fHitsCollection;
    G4double fTotalEnergyDeposit;
};

}

#endif
