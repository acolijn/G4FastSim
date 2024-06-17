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
/// \file B1/src/RunAction.cc
/// \brief Implementation of the B1::RunAction class

#include "RunAction.hh"
#include "EventAction.hh"	
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4AnalysisManager.hh"
#include "GammaRayHelper.hh"
#include <cmath>


namespace G4FastSim
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction* eventAction, GammaRayHelper* helper)
  : fEventAction(eventAction), fGammaRayHelper(helper)
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create the generic analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  InitializeNtuples();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::InitializeNtuples(){

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  if (G4Threading::IsMasterThread()) {
        analysisManager->SetDefaultFileType("root");
        analysisManager->OpenFile("G4FastSim");
        //  
        analysisManager->SetVerboseLevel(1);
        // Default settings
        analysisManager->SetNtupleMerging(true);
        // Creating event data ntuple
        DefineEventNtuple();
        // Creating and filling physics data ntuple
        DefineCrossSectionNtuple();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::DefineEventNtuple(){
  // Creating ntuple
  auto analysisManager = G4AnalysisManager::Instance();

  G4cout << "RunAction::BeginOfRunAction: Creating event data ntuple" << G4endl;

  eventNtupleId = analysisManager->CreateNtuple("ev", "G4FastSim ntuple");
  analysisManager->CreateNtupleDColumn(eventNtupleId, "Edep");     // column Id = 0
  analysisManager                                   // column Id = 1
    ->CreateNtupleDColumn(eventNtupleId, "xh", fEventAction->GetX());
  analysisManager                                   // column Id = 2
    ->CreateNtupleDColumn(eventNtupleId, "yh", fEventAction->GetY());
  analysisManager->FinishNtuple(eventNtupleId);
  G4cout <<"RunAction::BeginOfRunAction: Event data ntuple created. ID = "<< eventNtupleId << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::DefineCrossSectionNtuple(){
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  // Create ntuple for cross-section data
  G4cout << "RunAction::BeginOfRunAction: Creating cross-section data ntuple" << G4endl;

  crossSectionNtupleId = analysisManager->CreateNtuple("gam", "Gamma-ray cross-section data");
  analysisManager->CreateNtupleSColumn(crossSectionNtupleId, "mat"); // material
  analysisManager->CreateNtupleSColumn(crossSectionNtupleId, "proc"); // process
  analysisManager->CreateNtupleDColumn(crossSectionNtupleId, "e"); // energy
  analysisManager->CreateNtupleDColumn(crossSectionNtupleId, "att"); // attentuation length
  analysisManager->FinishNtuple(crossSectionNtupleId);
  G4cout <<"RunAction::BeginOfRunAction: Cross-section data ntuple created. ID = "<< crossSectionNtupleId << G4endl;

  //analysisManager->OpenFile();

  // Calculate the cross-sections and fill the HDF5 ntuple
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  G4Material* material = (*materialTable)[2];
  
  double startEnergy = 0.001 * MeV;
  double endEnergy = 10.0 * MeV;
  int numSteps = 1000;
  double factor = std::pow(endEnergy / startEnergy, 1.0 / (numSteps - 1));

  std::vector<G4String> processNames = {"compton", "phot", "tot", "att"};

  for (auto* mat : *materialTable) {
    for (auto processName : processNames) {
      for (int i = 0; i < numSteps; i++) {
        double energy = startEnergy * std::pow(factor, i);
        double crossSection = 0;
        if (processName == "compton") {
          crossSection = fGammaRayHelper->GetComptonCrossSection(energy, mat);
        } else if (processName == "phot") {
          crossSection = fGammaRayHelper->GetPhotoelectricCrossSection(energy, mat);
        } else if (processName == "tot") {
          crossSection = fGammaRayHelper->GetTotalCrossSection(energy, mat);
        }
        analysisManager->FillNtupleSColumn(crossSectionNtupleId, 0, mat->GetName());
        analysisManager->FillNtupleSColumn(crossSectionNtupleId, 1, processName);
        analysisManager->FillNtupleDColumn(crossSectionNtupleId, 2, energy / MeV);
        if (processName == "att") {
          analysisManager->FillNtupleDColumn(crossSectionNtupleId, 3, fGammaRayHelper->GetMassAttenuationCoefficient(energy, mat)/cm2);
        } else {
          analysisManager->FillNtupleDColumn(crossSectionNtupleId, 3, crossSection);
        }

        analysisManager->AddNtupleRow(crossSectionNtupleId);
      }
    }
    G4cout << mat->GetName() <<" density = "<< mat->GetDensity() / (g/cm3) <<" attenutation at 1 MeV: " << fGammaRayHelper->GetMassAttenuationCoefficient(1.0 * MeV, mat) / (cm2/g)  << " " << G4endl;
    G4double thickness = 1.0 * cm;
    G4double att = fGammaRayHelper->GetMassAttenuationCoefficient(1.0 * MeV, mat) * mat->GetDensity() * thickness;
    G4cout << mat->GetName() <<" linear attenuation at 1 MeV for 1 cm thickness: " << att << " " << G4endl;
  }
  G4cout << "units.... cm="<< cm << " MeV=" << MeV << " g= "<<g<<G4endl; 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{

  // save histograms & ntuple
  //
  if (G4Threading::IsMasterThread()) {
      auto analysisManager = G4AnalysisManager::Instance();
      analysisManager->Write();
      analysisManager->CloseFile();
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
