#include "SensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "Hit.hh"

namespace G4FastSim {

SensitiveDetector::SensitiveDetector(const G4String& name, const G4String& hitsCollectionName)
    : G4VSensitiveDetector(name), fHitsCollection(nullptr), fHitsCollectionID(-1), fTotalEnergyDeposit(0.) {
    collectionName.insert(hitsCollectionName);
}

SensitiveDetector::~SensitiveDetector() {}

void SensitiveDetector::Initialize(G4HCofThisEvent* hce) {
    
    fHitsCollection = new HitsCollection(SensitiveDetectorName, collectionName[0]);
    if (fHitsCollectionID < 0) {
        fHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
    }
    //fHitsCollection = new G4THitsCollection<Hit>(SensitiveDetectorName, collectionName[0]);

    G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(hcID, fHitsCollection);
    fTotalEnergyDeposit = 0.;  // Reset total energy deposit

}

G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {
    G4double edep = step->GetTotalEnergyDeposit();
    if (edep == 0.) return false;

    G4FastSim::Hit* newHit = new G4FastSim::Hit();
    //Hit* newHit = new Hit();
    newHit->energyDeposit = edep;
    newHit->position = step->GetPreStepPoint()->GetPosition();
    newHit->time = step->GetPreStepPoint()->GetGlobalTime();
    newHit->trackID = step->GetTrack()->GetTrackID();
    newHit->parentID = step->GetTrack()->GetParentID();
    newHit->momentum = step->GetPreStepPoint()->GetMomentum();
    newHit->particleType = step->GetTrack()->GetDefinition()->GetParticleName();

    fHitsCollection->insert(newHit);
    fTotalEnergyDeposit += edep;  // Accumulate total energy deposit

    return true;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent* hce) {
    if (!fHitsCollection) return;

    // Example of processing hits at the end of the event
    G4int numHits = fHitsCollection->entries();
    //G4cout << "End of event: total energy deposited in sensitive detector: "
    //       << fTotalEnergyDeposit / CLHEP::MeV << " MeV" << G4endl;
}

}
