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

  if (verbosityLevel >= 2){
    G4cout <<"Entering SteppingAction::UserSteppingAction"<<G4endl;
  }

  // get the energy of the particle
  G4double energy = step->GetPreStepPoint()->GetKineticEnergy();
  // get volume of the current step
  G4LogicalVolume* volume_pre = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4Material* material        = volume_pre->GetMaterial();
  // get the attenuation length of the material for the gamma ray at the current energy
  G4double attenuation_length = fGammaRayHelper->GetAttenuationLength(energy, material);
  G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

  // kill track if it is in the big bad world.....
  if(volume_pre->GetName() == "World") {
    step->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }

  if (verbosityLevel >= 2) Print(step);
  // if (1) we deal with the original gamma ray and (2) we are inside the fiducial volume and (3) we have not reached the maximum number of scatters ==> go scatter dude!
  G4int number_of_scatters = fEventAction->GetNumberOfScatters();
  if( (particleName == "geantino") &&
      (volume_pre->GetName() == "LXeFiducial") && 
      (number_of_scatters < fEventAction->GetNumberOfScattersMax())) {
    //
    // * Generate an interaction somewhere in the LXeFiducial volume
    //
    auto result = fGammaRayHelper->GenerateInteractionPoint(step);
    G4ThreeVector interactionPoint = result.first;
    // calculate the weight of the event and add it to teh event sum of logs
    fEventAction->AddWeight(std::log(result.second));

    // * scattering:
    //     - if the maximum allowed energy is above the photo-peak, generate a PE or Compton scatter, based on the relative cross sections
    //     - if the maximum allowed energy is below the photo-peak ->
    //                  i) ignore PE effect and assign a weight. Then just do Compton scatter
    //                  ii) calculate the maximum scattering angle possible and generate a Compton scatter. Calculate the event weight based on the non-sampled scattering angles
    G4double weight = DoScatter(step, interactionPoint);
    fEventAction->AddWeight(std::log(weight));

    // update the number of scatters
    fEventAction->SetNumberOfScatters(number_of_scatters + 1);

  } else {
    //  Update the event weight based on the traversed material. Once the geantino leaves the
    //  active vlume of the TPC the track will be killed and no further weights need to be calculated
    //  NOTE: at a later stage we could generate a gamma ray with the appropriate energy at the boundary
    //  of the active volume. This can facilitate generation of events that bounce back into the fiducial 
    //  volume
    
    // check if the particle is a geantino
    if (particleName == "geantino") {
      G4double weight = std::exp(-step->GetStepLength() / attenuation_length);
      fEventAction->AddWeight(std::log(weight));

      // if the next volume is the InnerCryostat and the geantino has undergone one or multiple scatters, 
      // we want to kill the track
      G4String nextVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
      if ((nextVolume == "InnerCryostat") && (number_of_scatters > 0)){
        step->GetTrack()->SetTrackStatus(fStopAndKill);
      }

    } else {
      // not a geantino.... so nothing special needs to be done here.
      return;
    }
    
  }


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/**
 * Performs scattering of particles based on their energy and material properties.
 * 
 * @param step The G4Step object containing information about the current step.
 * @param x0 The position vector of the current step.
 */
G4double SteppingAction::DoScatter(const G4Step* step, G4ThreeVector x0){
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // get the energy of the particle
  G4double energy = step->GetPreStepPoint()->GetKineticEnergy();
  G4double energyDeposit = 0.0;
  G4double weight = 1.0;

  // get the volume of the current step
  G4LogicalVolume* volume_pre = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4Material* material        = volume_pre->GetMaterial();

  G4double photoelectricCrossSection = fGammaRayHelper->GetPhotoelectricCrossSection(energy, material);
  G4double comptonCrossSection       = fGammaRayHelper->GetComptonCrossSection(energy, material);
  G4double ratio = comptonCrossSection / (photoelectricCrossSection + comptonCrossSection);

  G4double maxEnergy = fEventAction->GetAvailableEnergy();

  G4double rand = G4UniformRand();

  //G4cout << "DoScatter:: PE cross section: " << photoelectricCrossSection << G4endl;
  //G4cout << "DoScatter:: Compton cross section: " << comptonCrossSection << G4endl;
  //G4cout << "DoScatter:: Random number: " << rand << G4endl;
  
  if ( (rand > ratio) && (energy < maxEnergy) ){
    // generate a PE scatter: we make an energy deposit with the full energyand kill the track
    energyDeposit = energy;
  } else {
    // take into account that we ignored the PE effect and assign a weight
    // generate a Compton scatter and update the energy deposit
    InteractionData compton = fGammaRayHelper->DoComptonScatter(step, x0, maxEnergy);
    energyDeposit = compton.energyDeposited;

    if(energy > maxEnergy) {
      weight *= ratio; // take into account ignored PE effect
      weight *= compton.weight; // take into account the Compton scatter weight, only if the maximum allowed energy deposit is somewhere in the Compton region
    } 

    //
    // make a secondary track
    //
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDefinition = particleTable->FindParticle("geantino");
    // Create the dynamic particle
    G4DynamicParticle* dynamicParticle = new G4DynamicParticle(particleDefinition, compton.dir, compton.energy);
    G4Track* newTrack = new G4Track(dynamicParticle, 0.0, x0);
    // Set additional properties of the new track if needed

    G4Track* track = step->GetTrack();
    newTrack->SetParentID(track->GetTrackID());
    newTrack->SetGoodForTrackingFlag(true);
    // add new track to collection of secondaries
    G4TrackVector* secondaries = const_cast<G4TrackVector*>(step->GetSecondary());
    secondaries->push_back(newTrack);

    analysisManager->FillH1(0, compton.cosTheta);
  }

  // update the available energy
  fEventAction->ReduceAvailableEnergy(energyDeposit);  
  // kill the track!
  G4Track* track = step->GetTrack();
  track->SetTrackStatus(fStopAndKill);

    // create a new hit
  Hit* newHit = new Hit();
  newHit->energyDeposit = energyDeposit;
  newHit->position = x0;
  newHit->time = step->GetPreStepPoint()->GetGlobalTime();
  newHit->trackID = step->GetTrack()->GetTrackID();
  newHit->parentID = step->GetTrack()->GetParentID();
  newHit->particleType = "manual";

  AddHitToCollection(newHit, "LXeFiducialCollection");

  return weight;
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

  G4cout << "SteppingAction::Print track ID             : " << track->GetTrackID() << G4endl;
  G4cout << "                      track parent ID      : " << track->GetParentID() << G4endl;
  G4cout << "                      particle name        : " << track->GetParticleDefinition()->GetParticleName() << G4endl;
  G4cout << "                      energy               : " << track->GetKineticEnergy() / keV << " keV"<< G4endl;
  G4cout << "                      volume_pre name      : " << volume_pre->GetName() << G4endl;
  G4cout << "                      volume_pre position  : " << step->GetPreStepPoint()->GetPosition()/cm << " cm"<<G4endl;
  G4cout << "                      volume_post position : " << step->GetPostStepPoint()->GetPosition()/cm << " cm"<<G4endl;
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
        if (verbosityLevel >=2) G4cout << "Initializing hits collections" << G4endl;
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        fHitsCollectionIDs["LXeCollection"] = SDman->GetCollectionID("LXeCollection");
        fHitsCollectionIDs["LXeFiducialCollection"] = SDman->GetCollectionID("LXeFiducialCollection");
        if (verbosityLevel >=2){
          G4cout << "Hits collection ID for LXeCollection: " << fHitsCollectionIDs["LXeCollection"] << G4endl;
          G4cout << "Hits collection ID for LXeFiducialCollection: " << fHitsCollectionIDs["LXeFiducialCollection"] << G4endl;
        }
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
