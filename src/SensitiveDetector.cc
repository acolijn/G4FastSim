#include "SensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4VProcess.hh"
#include "Hit.hh"

namespace G4Sim {

SensitiveDetector::SensitiveDetector(const G4String& name, const G4String& volumeName)
    : G4VSensitiveDetector(name), fTotalEnergyDeposit(0.) {
    
    collectionName.insert(volumeName + "Collection");  // Register hit collection for each volume
    
}

void SensitiveDetector::Initialize(G4HCofThisEvent* hce) {
    fHitsCollections.clear();  // Clear any previous hit collections

    // Initialize a hit collection for each registered name
    for (size_t i = 0; i < collectionName.size(); ++i) {
        G4String collectionNameStr = collectionName[i];
        auto* hitsCollection = new HitsCollection(SensitiveDetectorName, collectionNameStr);
        fHitsCollections.push_back(hitsCollection);

        G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionNameStr);
        if (hcID < 0) {
            hcID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
        }

        hce->AddHitsCollection(hcID, hitsCollection);  // Register each collection with the event
    }

    fTotalEnergyDeposit = 0.0;  // Reset total energy deposit
}

// void SensitiveDetector::Initialize(G4HCofThisEvent* hce) {
    
//     fHitsCollection = new HitsCollection(SensitiveDetectorName, collectionName[0]);
//     if (fHitsCollectionID < 0) {
//         fHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
//     }
//     //fHitsCollection = new G4THitsCollection<Hit>(SensitiveDetectorName, collectionName[0]);

//     G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
//     hce->AddHitsCollection(hcID, fHitsCollection);
//     fTotalEnergyDeposit = 0.;  // Reset total energy deposit

// }

SensitiveDetector::~SensitiveDetector() {}


G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {
    G4double edep = step->GetTotalEnergyDeposit();
    if (edep == 0.) return false;

    // Get the volume where the step occurred
    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    //G4TouchableHandle touchable = preStepPoint->GetTouchable();
    const G4VTouchable* touchable = preStepPoint->GetTouchable();  // Keep it const since we are not modifying it

    G4String volumeName = touchable->GetVolume()->GetName();  // Get the name of the current volume

    // Loop through the registered hit collections and find the correct one for the volume
    for (size_t i = 0; i < fHitsCollections.size(); ++i) {
        G4String name = volumeName + "Collection";
        if (name == collectionName[i]) {  // Match the hit collection with the volume
            auto* newHit = new G4Sim::Hit();
            newHit->energyDeposit = edep;
            newHit->position = step->GetPostStepPoint()->GetPosition();
            newHit->time = step->GetPostStepPoint()->GetGlobalTime();
            newHit->trackID = step->GetTrack()->GetTrackID();
            newHit->parentID = step->GetTrack()->GetParentID();
            newHit->momentum = step->GetPreStepPoint()->GetMomentum();
            newHit->particleType = step->GetTrack()->GetDefinition()->GetParticleName();
            newHit->processType = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
            newHit->particleEnergy0 = step->GetPreStepPoint()->GetKineticEnergy();
            newHit->particleEnergy1 = step->GetPostStepPoint()->GetKineticEnergy();

            // Insert the hit into the correct hit collection
            fHitsCollections[i]->insert(newHit);
            break;  // We found the right collection, no need to continue the loop
        }
    }

    fTotalEnergyDeposit += edep;  // Accumulate total energy deposit
    return true;
}

void SensitiveDetector::AddCollection(G4String volumeName) {
    collectionName.insert(volumeName + "Collection");  // Register hit collection for each volume
}


void SensitiveDetector::EndOfEvent(G4HCofThisEvent* hce) {
    //if (!fHitsCollection) return;

    // Example of processing hits at the end of the event
    //G4int numHits = fHitsCollection->entries();
    //G4cout << "End of event: total energy deposited in sensitive detector: "
    //       << fTotalEnergyDeposit / CLHEP::MeV << " MeV" << G4endl;
}

} // namespace G4Sim
