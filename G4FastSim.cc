#include "DetectorConstruction.hh"
#include "DetectorConstructionMessenger.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "QBBC.hh"
#include "FTFP_BERT_HP.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4PhysListFactory.hh"
#include "G4EmLivermorePhysics.hh"

#include "Randomize.hh"
#include "GammaRayHelper.hh"
#include "CustomEmPhysics.hh"

using namespace G4FastSim;

/**
 * @brief The main function of the program.
 *
 * This function is the entry point of the program. It initializes the necessary components,
 * sets up the run manager, initializes the visualization, and processes the macro or starts
 * the UI session based on the command line arguments. After the execution, it frees the memory
 * allocated for the visualization manager and the run manager.
 *
 * @param argc The number of command line arguments.
 * @param argv An array of command line arguments.
 * @return An integer representing the exit status of the program.
 */
int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);

  //use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // Construct the default run manager
  //
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  //runManager->SetNumberOfThreads(1);

  // Set mandatory initialization classes
  //
  // Detector construction

  GammaRayHelper* helper = &GammaRayHelper::Instance();
  runManager->SetUserInitialization(new DetectorConstruction(helper));

  // Physics list
  G4PhysListFactory factory;
  G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("FTFP_BERT_HP");
  physicsList->ReplacePhysics(new G4EmLivermorePhysics());
  //if you want to mess with the Em physics list ...... physicsList->ReplacePhysics(new CustomEmPhysics());
  runManager->SetUserInitialization(physicsList);
  // User action initialization
  runManager->SetUserInitialization(new ActionInitialization(helper));

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....