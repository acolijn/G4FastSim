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
/// \file B1/include/EventAction.hh
/// \brief Definition of the B1::EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>
#include <mutex>

/// Event action class
///

namespace G4FastSim
{

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    std::vector<G4double>& GetX(){return fX;};
    std::vector<G4double>& GetY(){return fY;};

    void AddEdep(G4double edep) { fEdep += edep; }

  private:
    // define here all the variables that you want to store for each event in the 
    // ntuple tree  
    static thread_local G4double fEdep;

    static thread_local std::vector<G4double> fEd;
    static thread_local std::vector<G4double> fX;
    static thread_local std::vector<G4double> fY;

    static std::mutex mtx; // Mutex for thread safety

};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


