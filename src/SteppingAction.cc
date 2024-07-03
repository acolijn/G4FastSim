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
#include "GammaRayHelper.hh"

#include <vector>
#include <mutex>

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
      fGammaRayHelper(helper) {
      
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
  //if (!fScoringVolume) {
  //  const auto detConstruction = static_cast<const DetectorConstruction*>(
  //    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  //  fScoringVolume = detConstruction->GetScoringVolume();
  //}

  G4Track* track = step->GetTrack();

  // get volume of the current step
  G4LogicalVolume* volume_pre = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4Material* material = volume_pre->GetMaterial();

  auto poststep = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
  // if you leave the Geant4 World volume.... do nothing
  G4LogicalVolume* volume_post = nullptr;
  if(poststep) {
    volume_post  = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  }

  // if Verbose level >=2 print this information
  if(verbosityLevel >= 2 )
  {
    G4cout << "track ID: " << track->GetTrackID() << G4endl;
    G4cout << "  volume_pre: " << volume_pre->GetName() << G4endl;
    G4cout << "  material: " << material->GetName() << " Attenuation length:" <<  fGammaRayHelper->GetAttenuationLength(1.0 * MeV, material) /cm << " cm" << G4endl;
  }
  
  // get the track information

  G4ParticleDefinition* particle = track->GetDefinition();
  G4String particleName = particle->GetParticleName();

  // change direction of teh track
  //G4ThreeVector direction = track->GetMomentumDirection();
  //G4ThreeVector newDirection = direction;
  //newDirection.setY(direction.getY()+0.1);
  //newDirection = newDirection.unit();
  //track->SetMomentumDirection(newDirection);

  //G4cout<< "Track direction: " << direction << G4endl;
  //G4cout<< "New track direction: " << newDirection << G4endl;

  // get the particle information

  //G4cout << "Particle name: " << particleName << G4endl;

  // set the particle definition to gamma

  //G4ParticleDefinition* particle_new = particleTable->FindParticle("gamma");
  //if(volume_post->GetName() == "LXe") {
  //}
  //  track->SetDefinition(particle_new);
  //  //track->SetParticleEnergy(3.*MeV);
  //}

  //G4cout << "Step length: " << step->GetStepLength() << G4endl;

//    step->AddSecondaryParticle(newTrack);

    //AddSecondaryParticle(newTrack);
 
  if((step->GetTrack()->GetTrackID() == 1) && (volume_pre->GetName() == "LXe")) {
    // Get the pre-step point
    G4StepPoint* preStepPoint = step->GetPreStepPoint();

    // Get the position and time of the pre-step point
    G4ThreeVector position = preStepPoint->GetPosition();
    G4double globalTime = preStepPoint->GetGlobalTime();

    // Get the particle table
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

    // Define the particle type (e.g., proton)
    G4ParticleDefinition* particleDefinition = particleTable->FindParticle("geantino");

    if (particleDefinition == nullptr) {
        G4cerr << "Error: Particle type not found." << G4endl;
        return;
    }

    // Define the momentum direction and kinetic energy of the new particle
    G4ThreeVector momentumDirection(1.0, 0.0, 1.0); // example direction
    momentumDirection = momentumDirection.unit();
    G4double kineticEnergy = 5.0 * MeV; // example energy

    // Create the dynamic particle
    G4DynamicParticle* dynamicParticle = new G4DynamicParticle(particleDefinition, momentumDirection, kineticEnergy);

    // Create the new track
    G4Track* newTrack = new G4Track(dynamicParticle, globalTime, position);

    // Set additional properties of the new track if needed
    newTrack->SetParentID(step->GetTrack()->GetTrackID());
    newTrack->SetGoodForTrackingFlag(true);

    // Add the new track to the list of secondaries
    G4TrackVector* secondaries = const_cast<G4TrackVector*>(step->GetSecondary());

    secondaries->push_back(newTrack);

    track->SetTrackStatus(fStopAndKill);
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
