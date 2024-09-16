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
#include "G4Gamma.hh"
#include "G4PhysicalConstants.hh"

#include <cmath>

namespace G4Sim{
/**
 * @file RunAction.cc
 * @brief Implementation of the RunAction class.
 *
 * The RunAction class is responsible for managing actions that occur during a run of the simulation.
 * It initializes and defines the analysis manager, creates and fills ntuples for event data, cross-section data,
 * and differential cross-section data, and performs actions at the beginning and end of a run.
 */


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction* eventAction, GammaRayHelper* helper)
  : fEventAction(eventAction), fGammaRayHelper(helper)
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1000);

  // Create the generic analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  fMessenger = new RunActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{

  // Get the initial energy from the primary generator action
  const PrimaryGeneratorAction* primaryGeneratorAction = static_cast<const PrimaryGeneratorAction*>(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

  G4cout << "Runaction::BeginOfRunAction: E0 = " << primaryGeneratorAction->GetInitialEnergy() / keV << " keV" << G4endl;
  // Initialize the gamma-ray helper
  fGammaRayHelper->Initialize();

  // initialize the analysis manager and ntuples
  InitializeNtuples();

  // passing some parameters to the event action
  G4cout <<"RunAction::BeginOfRunAction: Normal (0) - or - FastSimulation (1) = "<< fFastSimulation << G4endl;
  fEventAction->SetFastSimulation(fFastSimulation);
  G4cout <<"RunAction::BeginOfRunAction: Maximum number of scatters = "<< fNumberOfScattersMax << G4endl;
  fEventAction->SetNumberOfScattersMax(fNumberOfScattersMax);
  G4cout <<"RunAction::BeginOfRunAction: Maximum energy deposit = "<< fMaxEnergy << G4endl;
  fEventAction->SetMaxEnergy(fMaxEnergy);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::~RunAction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::InitializeNtuples(){

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  if (G4Threading::IsMasterThread()) {
    analysisManager->SetDefaultFileType("root");
    analysisManager->OpenFile(fOutputFileName);
    //  
    analysisManager->SetVerboseLevel(1);
    // Default settings
    analysisManager->SetNtupleMerging(true);

    // Creating histograms
    analysisManager->CreateH1("cost", "cos theta of Compton", 2200, -1.1, +1.1); // id = 0

    // Creating event data ntuple
    DefineEventNtuple();
    // Creating and filling physics data ntuple
    DefineCrossSectionNtuple();
    // Creating and filling differential cross-section data ntuple
    // done from EventAction at the first event: since we the know the energy of the gamma rays. DefineDifferentialCrossSectionNtuple();
  }
}

/**
 * @brief Defines the differential cross-section Ntuple.
 * 
 * This function creates an Ntuple to store the differential cross-section data. It retrieves the material table,
 * creates the necessary columns in the Ntuple, and fills the Ntuple with the differential cross-section values
 * for each element in each material. The Ntuple columns include the material name, scattering angle, form factor,
 * Klein-Nishina cross-section, and atomic number of the element.
 */
void RunAction::DefineDifferentialCrossSectionNtuple(G4double e0) const {
  auto analysisManager = G4AnalysisManager::Instance();

  // get material table
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();

  diffXsecNtupleId = analysisManager->CreateNtuple("diff_xsec", "differential cross-section data");
  analysisManager->CreateNtupleSColumn(diffXsecNtupleId, "mat");    // column Id = 0
  analysisManager->CreateNtupleDColumn(diffXsecNtupleId, "cost");     // column Id = 1
  analysisManager->CreateNtupleDColumn(diffXsecNtupleId, "ff");     // column Id = 2
  analysisManager->CreateNtupleDColumn(diffXsecNtupleId, "kn");     // column Id = 3
  analysisManager->CreateNtupleIColumn(diffXsecNtupleId, "Z");     // column Id = 4
  analysisManager->FinishNtuple(diffXsecNtupleId);

  G4cout <<"RunAction::BeginOfRunAction: Diff Xsec ntuple created. ID = "<< diffXsecNtupleId << G4endl;
  G4cout <<"RunAction::BeginOfRunAction: Material table size = "<< materialTable->size() << G4endl;
  G4cout <<"RunAction::BeginOfRunAction: create ntuple for energy = " << e0/MeV << " MeV" << G4endl;


  std::vector<std::string> elements_used;
  for (size_t i = 0; i < materialTable->size(); ++i) {

    G4Material* material = (*materialTable)[i];
    //auto* comptonModel = fGammaRayHelper->GetComptonModel(material);
    G4MaterialCutsCouple* cuts = new G4MaterialCutsCouple(material, 0);

    G4DynamicParticle* gamma = new G4DynamicParticle(G4Gamma::Gamma(), G4ThreeVector(1,0,0), e0);
    G4ParticleDefinition *particle = gamma->GetDefinition();

    // sloop over all elements in a material
    const G4ElementVector* elementVector = material->GetElementVector();
    for (size_t i = 0; i < material->GetNumberOfElements(); ++i) {
      const G4Element* elm = (*elementVector)[i];
      if (std::find(elements_used.begin(), elements_used.end(), elm->GetName()) == elements_used.end()) {

        for (G4double cost = -1.0; cost <= 1.0; cost += 0.001) {
          // get scatter function
          G4double ff = fGammaRayHelper->GetComptonModel()->FormFactor(elm, gamma, cost);
          // get differential cross-section (Klein-Nishina)
          G4double kn = fGammaRayHelper->GetComptonModel()->KleinNishina(gamma, cost);

          analysisManager->FillNtupleSColumn(diffXsecNtupleId, 0, elm->GetName());
          analysisManager->FillNtupleDColumn(diffXsecNtupleId, 1, cost);
          analysisManager->FillNtupleDColumn(diffXsecNtupleId, 2, ff);
          analysisManager->FillNtupleDColumn(diffXsecNtupleId, 3, kn);
          analysisManager->FillNtupleIColumn(diffXsecNtupleId, 4, elm->GetZ());
          analysisManager->AddNtupleRow(diffXsecNtupleId);
        }
        elements_used.push_back(elm->GetName());
      }
    }
  }
}
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Defines the event ntuple for data analysis.
 * 
 * This function creates an event ntuple using the G4AnalysisManager class. The ntuple contains columns for storing energy deposition (Edep), x-coordinate (xh), and y-coordinate (yh) of each event. The ntuple is finished and assigned an ID.
 */
void RunAction::DefineEventNtuple(){
  // Creating ntuple
  auto analysisManager = G4AnalysisManager::Instance();

  G4cout << "RunAction::BeginOfRunAction: Creating event data ntuple" << G4endl;

  eventNtupleId = analysisManager->CreateNtuple("ev", "G4FastSim ntuple");
  analysisManager->CreateNtupleDColumn(eventNtupleId, "ev");   // column Id = 0
  analysisManager->CreateNtupleDColumn(eventNtupleId, "w");    // column Id = 1
  analysisManager->CreateNtupleDColumn(eventNtupleId, "type"); // column Id = 2
  analysisManager->CreateNtupleDColumn(eventNtupleId, "xp");   // column Id = 3
  analysisManager->CreateNtupleDColumn(eventNtupleId, "yp");   // column Id = 4
  analysisManager->CreateNtupleDColumn(eventNtupleId, "zp");   // column Id = 5
  analysisManager->CreateNtupleDColumn(eventNtupleId, "eh", fEventAction->GetE()); 
  analysisManager->CreateNtupleDColumn(eventNtupleId, "xh", fEventAction->GetX()); 
  analysisManager->CreateNtupleDColumn(eventNtupleId, "yh", fEventAction->GetY()); 
  analysisManager->CreateNtupleDColumn(eventNtupleId, "zh", fEventAction->GetZ()); 
  analysisManager->CreateNtupleDColumn(eventNtupleId, "wh", fEventAction->GetW()); 
  analysisManager->CreateNtupleIColumn(eventNtupleId, "id", fEventAction->GetID()); 
  analysisManager->CreateNtupleDColumn(eventNtupleId, "edet", fEventAction->GetEdet());
  analysisManager->CreateNtupleIColumn(eventNtupleId, "ndet", fEventAction->GetNdet());
  analysisManager->CreateNtupleIColumn(eventNtupleId, "nphot", fEventAction->GetNphot());
  analysisManager->CreateNtupleIColumn(eventNtupleId, "ncomp", fEventAction->GetNcomp());
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
  
  double startEnergy =  1 * keV;
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
          analysisManager->FillNtupleDColumn(crossSectionNtupleId, 3, crossSection / barn);
        }

        analysisManager->AddNtupleRow(crossSectionNtupleId);
      }
    }
    //G4cout << mat->GetName() <<" density = "<< mat->GetDensity() / (g/cm3) <<" attenutation at 1 MeV: " << fGammaRayHelper->GetMassAttenuationCoefficient(1.0 * MeV, mat) / (cm2/g)  << " " << G4endl;
    //G4double thickness = 1.0 * cm;
    //G4double att = fGammaRayHelper->GetMassAttenuationCoefficient(1.0 * MeV, mat) * mat->GetDensity() * thickness;
    //G4cout << mat->GetName() <<" linear attenuation at 1 MeV for 1 cm thickness: " << att << " " << G4endl;
  }
  //G4cout << "units.... cm="<< cm << " MeV=" << MeV << " g= "<<g<<G4endl; 
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
