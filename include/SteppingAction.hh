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
/// \file B1/include/SteppingAction.hh
/// \brief Definition of the B1::SteppingAction class

#ifndef B1SteppingAction_h
#define B1SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"
#include "GammaRayHelper.hh"
#include "Hit.hh"
#include "SensitiveDetector.hh"
#include <map>
#include <vector>
#include <utility>

class G4LogicalVolume;

/// Stepping action class
///

namespace G4FastSim
{

class EventAction;

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(EventAction* eventAction, GammaRayHelper* helper);
    ~SteppingAction();// override = default;

    // method from the base class
    void UserSteppingAction(const G4Step*) override;
    void SetVerbosity(G4int level) { verbosityLevel = level; }  
    void Print(const G4Step* step);

  private:
    void AddHitToCollection(Hit *newhit, G4String collectionName);
    G4double DoScatter(const G4Step* step, G4ThreeVector x0);


    EventAction* fEventAction;
    G4LogicalVolume* fScoringVolume;
    G4ParticleTable* particleTable;
    GammaRayHelper* fGammaRayHelper;

    std::map<G4String, HitsCollection*> fHitsCollections;
    std::map<G4String, G4int> fHitsCollectionIDs;
    G4bool fHitsCollectionInitialized;

    G4int verbosityLevel=0;

};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
