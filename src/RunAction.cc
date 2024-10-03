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

/**
 * @namespace G4Sim
 * @brief Namespace for the G4Sim library.
/*/
namespace G4Sim{
/**
 * @file RunAction.cc
 * @brief Implementation of the RunAction class.
 *
 * The RunAction class is responsible for managing actions that occur during a run of the simulation.
 * It initializes and defines the analysis manager, creates and fills ntuples for event data, cross-section data,
 * and differential cross-section data, and performs actions at the beginning and end of a run.
 */
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

/**
 * @brief This method is called at the beginning of each run.
 *
 * It performs the following actions:
 * - Retrieves the initial energy from the primary generator action and prints it.
 * - Initializes the gamma-ray helper.
 * - Initializes the analysis manager and ntuples.
 * - Passes parameters to the event action, including:
 *   - Fast simulation mode (normal or fast).
 *   - Maximum number of scatters.
 *   - Maximum energy deposit.
 *
 * @param run Pointer to the current G4Run object.
 */
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

/**
 * @brief Initializes the ntuples for data analysis.
 *
 * This function sets up the analysis manager, configures the output file type,
 * opens the output file, and sets various default settings. It also creates
 * histograms and defines ntuples for event data, physics data, and differential
 * cross-section data.
 *
 * The function performs the following steps:
 * - Retrieves the analysis manager instance.
 * - If the current thread is the master thread:
 *   - Sets the default file type to "root".
 *   - Opens the output file with the specified name.
 *   - Sets the verbose level to 1.
 *   - Enables ntuple merging.
 *   - Creates a histogram for the cosine of the Compton scattering angle.
 *   - Defines the ntuple for event data.
 *   - Defines the ntuple for physics data.
 *   - Prepares for the definition of the differential cross-section ntuple,
 *     which will be done in the EventAction at the first event.
 */
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
 * This function creates an event ntuple using the G4AnalysisManager class. 
 * 
 * Definition of the ntuple columns:
 * - ev = Event number
 * - w  = Event weight (only for fast simulation)
 * - type = Particle type (only for fast simulation: did the gamma scatter outside the detector prior to makingan interaction in the xenon?)
 * - xp = x position of the primary event
 * - yp = y position of the primary event
 * - zp = z position of the primary event
 * - eh = vector of energy deposited in the detector
 * - xh = vector of x position of the clusters in the detector
 * - yh = vector of y position of the clusters in the detector
 * - zh = vector of z position of the clusters in the detector
 * - wh = vector of weights (all the same, only for fast simulation)
 * - id = vector of detector IDs of the clusters
 * - edet = total energy deposited in each detector
 * - ndet = number of clusters in each detector
 * - nphot = number of photo-electric interactions in each detector
 * - ncomp = number of Compton interactions in each detector
 */
void RunAction::DefineEventNtuple(){
  // Creating ntuple
  auto analysisManager = G4AnalysisManager::Instance();

  G4cout << "RunAction::BeginOfRunAction: Creating event data ntuple" << G4endl;

  eventNtupleId = analysisManager->CreateNtuple("ev", "G4Sim ntuple");
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

/**
 * @brief Defines and fills a ntuple with gamma-ray cross-section data for various materials and processes.
 *
 * This function initializes an ntuple to store gamma-ray cross-section data, including material names,
 * process names, energy levels, and cross-section values. It iterates over a predefined set of materials
 * and processes, calculates the cross-sections for a range of energy levels, and fills the ntuple with
 * the computed data.
 *
 * The ntuple contains the following columns:
 * - "mat": Material name (string)
 * - "proc": Process name (string)
 * - "e": Energy (double, in MeV)
 * - "att": Attenuation length or cross-section (double, in cm^2/g or barns)
 *
 * The processes considered are:
 * - "compton": Compton scattering
 * - "phot": Photoelectric effect
 * - "tot": Total cross-section
 * - "att": Attenuation coefficient
 *
 * The energy range for the calculations is from 1 keV to 10 MeV, divided into 1000 logarithmic steps.
 *
 * @note This function assumes that the G4AnalysisManager and fGammaRayHelper instances are properly initialized.
 */
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
          analysisManager->FillNtupleDColumn(crossSectionNtupleId, 3, fGammaRayHelper->GetAttenuationLength(energy, mat)/cm);
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

/**
 * @brief Called at the end of a run to perform final actions.
 *
 * This method is executed at the end of a run. If the current thread is the master thread,
 * it saves histograms and ntuple data by writing and closing the analysis manager's file.
 *
 * @param run Pointer to the G4Run object representing the current run.
 */
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
