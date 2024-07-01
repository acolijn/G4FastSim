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
/// \file B1/src/EventAction.cc
/// \brief Implementation of the G4FastSim::EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

#include "GammaRayHelper.hh"



namespace G4FastSim
{
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::mutex EventAction::mtx;


EventAction::EventAction()
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{
  G4cout << "EventAction::BeginOfEventAction..... NEXT" << G4endl;	
  fEdep = 1.2345;
  fX.clear();
  fY.clear();

  //G4cout<<"EventAction::BeginOfEventAction next event...."<<G4endl;
  G4PrimaryVertex* primaryVertex = event->GetPrimaryVertex();
  fXp = primaryVertex->GetPosition().x();
  fYp = primaryVertex->GetPosition().y();
  fZp = primaryVertex->GetPosition().z();



  auto def =  event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition();
  G4cout << " def = " << def->GetParticleName() << G4endl;
  G4cout << " p = "<<event->GetPrimaryVertex()->GetPrimary()->GetMomentumDirection() << G4endl;
  // Kill the event if the particle does not point to the fiducial volume
  //if (TBD) {
  //  G4cout<<"EventAction::BeginOfEventAction: Killing event with fZp = "<<fZp<<G4endl;
  //  G4RunManager::GetRunManager()->AbortEvent();
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  G4cout << "EventAction::EndOfEventAction..... " << G4endl;
  // accumulate statistics in run action
  fX.push_back(0.1);
  fY.push_back(0.2);

  fX.push_back(0.3);
  fY.push_back(0.4);
  
  fX.push_back(-1.1);
  fY.push_back(-1.2);


  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Fill histograms
  // Fill ntuple
  //G4cout<<"EventAction::EndOfEventAction: fEdep = "<<fEdep<<G4endl;

    // Protect the following section with a mutex to ensure thread safety
  //{
  //    G4AutoLock lock(&mtx); // Use G4AutoLock for thread-safe initialization
  //std::lock_guard<std::mutex> lock(mtx);
  analysisManager->FillNtupleDColumn(0, 0, fEdep);
  analysisManager->FillNtupleDColumn(0, 1, fXp);
  analysisManager->FillNtupleDColumn(0, 2, fYp);
  analysisManager->FillNtupleDColumn(0, 3, fZp);
  analysisManager->AddNtupleRow(0);

  //G4cout<<"EventAction::EndOfEventAction: fX.size() = "<<fX.size()<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
