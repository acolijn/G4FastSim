//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4AutoLock.hh"
#include "G4AnalysisManager.hh"
#include "G4SDManager.hh"
#include "GammaRayHelper.hh"

#include <vector>
#include <mutex>
#include <utility>

namespace {
    G4Mutex mutex = G4MUTEX_INITIALIZER;
}

namespace G4FastSim
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction, GammaRayHelper* helper)
    : G4UserSteppingAction(),
      fEventAction(eventAction),
      fScoringVolume(nullptr),
      particleTable(G4ParticleTable::GetParticleTable()),
      fGammaRayHelper(helper), 
      fHitsCollectionInitialized(false){
  // the hits collection is not yet initialized
  fHitsCollectionInitialized = false;

}


SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fEventAction->IsFastSimulation()) {
    return;
  }
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


  if (verbosityLevel >= 2){
    G4cout <<"Entering SteppingAction::UserSteppingAction"<<G4endl;
  }

  // get the energy of the particle
  G4double energy = step->GetPreStepPoint()->GetKineticEnergy();
  // get volume of the current step
  G4LogicalVolume* volume_pre = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4Material* material        = volume_pre->GetMaterial();
  G4double attenuation_length = fGammaRayHelper->GetAttenuationLength(energy, material);

  // get the entrance and exit points of the step in the current volume
  G4ThreeVector entrance = step->GetPreStepPoint()->GetPosition();
  G4ThreeVector exit = step->GetPostStepPoint()->GetPosition();

  // get the attenuation length of the material for the gamma ray at the current energy
  G4double globalTime = step->GetPreStepPoint()->GetGlobalTime();

  Print(step);

  // if (1) we deal with the original gamma ray and (2) we are inside the fiducial volume and (3) we have not reached the maximum number of scatters ==> go scatter dude!
  G4int number_of_scatters = fEventAction->GetNumberOfScatters();
  if( (step->GetTrack()->GetTrackID() == 1) && 
      (volume_pre->GetName() == "LXeFiducial") && 
      (number_of_scatters < fEventAction->GetNumberOfScattersMax())) {
    //
    // * Generate an interaction somewhere in the LXeFiducial volume
    //
    auto result = GenerateInteractionPoint(entrance, exit, attenuation_length);
    G4ThreeVector interactionPoint = result.first;
    G4double weight = result.second;

    // * scattering:
    //     - if the maximum allowed energy is above the photo-peak, generate a PE or Compton scatter, based on the relative cross sections
    //     - if the maximum allowed energy is below the photo-peak ->
    //                  i) ignore PE effect and assign a weight. Then just do Compton scatter
    //                  ii) calculate the maximum scattering angle possible and generate a Compton scatter. Calculate the event weight based on the non-sampled scattering angles



    // interaction point
    G4ThreeVector momentumDirection = G4ThreeVector(0., 0., 1.);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDefinition = particleTable->FindParticle("geantino");
    // Create the dynamic particle
    G4DynamicParticle* dynamicParticle = new G4DynamicParticle(particleDefinition, momentumDirection, energy);
    // Create the new track
    G4Track* newTrack = new G4Track(dynamicParticle, globalTime, interactionPoint);
    // Set additional properties of the new track if needed
    G4Track* track = step->GetTrack();
    newTrack->SetParentID(track->GetTrackID());
    newTrack->SetGoodForTrackingFlag(true);

    // Add the new track to the list of secondaries
    G4TrackVector* secondaries = const_cast<G4TrackVector*>(step->GetSecondary());
    secondaries->push_back(newTrack);
    track->SetTrackStatus(fStopAndKill);


    // Create a new hit based on the energy deposit in the scattering event

    // this is just for testing .......
    Hit* newHit = new Hit();
    newHit->energyDeposit = 10.* MeV;
    newHit->position = interactionPoint;
    newHit->time = 1. ;
    newHit->trackID = -1;
    newHit->parentID = -1;
    newHit->momentum = momentumDirection;
    newHit->particleType = "manual";

    // the event to the hits collection
    AddHitToCollection(newHit, "LXeFiducialCollection");
    // update the number of scatters
    fEventAction->SetNumberOfScatters(number_of_scatters + 1);
  } else {
    // 1. update the event weight based on the traversed material

    // testing...
    // // //fGammaRayHelper->GenerateComptonScatteringDirection(volume_pre->GetMaterial(), step);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * Generates an interaction point along a line segment between two points.
 *
 * @param entrance The entrance point of the line segment.
 * @param exit The exit point of the line segment.
 * @param attenuation_length The attenuation length of the material.
 * @return The generated interaction point.
 */
std::pair<G4ThreeVector, G4double> SteppingAction::GenerateInteractionPoint(G4ThreeVector entrance, G4ThreeVector exit, G4double attenuation_length) {

  // maximum possible distance....
  G4double maxDistance = (exit - entrance).mag();

  // generate a random distance along the line segment
  G4double rand = G4UniformRand();
  G4double maxCDF = 1 - std::exp(-maxDistance / attenuation_length);
  G4double adjustedRand = rand * maxCDF;
  G4double distance = -attenuation_length * std::log(1 - adjustedRand);
  G4ThreeVector interactionPoint = entrance + (exit - entrance).unit() * distance;

  return std::make_pair(interactionPoint,maxCDF);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * Prints information about the given step.
 *
 * @param step The step for which to print information.
 */
void SteppingAction::Print(const G4Step* step) {

  G4Track* track = step->GetTrack();
  G4LogicalVolume* volume_pre = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  // if Verbose level >=2 print this information

  G4cout << "track ID        : " << track->GetTrackID() << G4endl;
  G4cout << "track parent ID : " << track->GetParentID() << G4endl;
  G4cout << "  particle name        : " << track->GetParticleDefinition()->GetParticleName() << G4endl;
  G4cout << "  energy               : " << track->GetKineticEnergy() / keV << " keV"<< G4endl;
  G4cout << "  volume_pre name      : " << volume_pre->GetName() << G4endl;
  G4cout << "  volume_pre position  : " << step->GetPreStepPoint()->GetPosition()/cm << " cm"<<G4endl;
  G4cout << "  volume_post position : " << step->GetPostStepPoint()->GetPosition()/cm << " cm"<<G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Adds a hit to the specified hits collection.
 *
 * This function adds the given hit to the hits collection with the specified collection name.
 * If the hits collection has not been initialized, it initializes the hits collection IDs.
 * If the hits collection ID is not found for the given collection name, an error message is printed and the function returns.
 * If the hits collection ID is -1 for the given collection name, an error message is printed and the function returns.
 * If the hits collection is not found for the given collection name, it retrieves the hits collection from the current event.
 * If the hits collection is still not found, an error message is printed and the function returns.
 *
 * @param newHit The hit to be added to the hits collection.
 * @param collectionName The name of the hits collection.
 */
void SteppingAction::AddHitToCollection(Hit* newHit, G4String collectionName){

     if (!fHitsCollectionInitialized) {
        G4cout << "Initializing hits collections" << G4endl;
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        fHitsCollectionIDs["LXeCollection"] = SDman->GetCollectionID("LXeCollection");
        fHitsCollectionIDs["LXeFiducialCollection"] = SDman->GetCollectionID("LXeFiducialCollection");
        G4cout << "Hits collection ID for LXeCollection: " << fHitsCollectionIDs["LXeCollection"] << G4endl;
        G4cout << "Hits collection ID for LXeFiducialCollection: " << fHitsCollectionIDs["LXeFiducialCollection"] << G4endl;
        // Add more collections here if needed
        fHitsCollectionInitialized = true;
    }

    if (fHitsCollectionIDs.find(collectionName) == fHitsCollectionIDs.end()) {
        G4cerr << "Error: Hits collection ID not found for collection name: " << collectionName << G4endl;
        return;
    }

    G4int hcID = fHitsCollectionIDs[collectionName];

    if (hcID == -1) {
        G4cerr << "Error: Hits collection ID is -1 for collection name: " << collectionName << G4endl;
        return;
    }

    const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
    G4HCofThisEvent* HCE = evt->GetHCofThisEvent();

    if (HCE) {
        fHitsCollections[collectionName] = static_cast<G4THitsCollection<Hit>*>(HCE->GetHC(hcID));
    }

    if (fHitsCollections[collectionName] == nullptr) {
        G4cerr << "Error: Hits collection not found for collection name: " << collectionName << G4endl;
        return;
    }

    //
    // Add the hit to the collection
    //
    fHitsCollections[collectionName]->insert(newHit);

}


}
