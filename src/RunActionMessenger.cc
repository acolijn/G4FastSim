#include "RunActionMessenger.hh"
#include "RunAction.hh"

///using namespace G4Sim;
namespace G4Sim {

RunActionMessenger::RunActionMessenger(RunAction* action)
    : fRunAction(action) {

    fFastSimulationCmd = new G4UIcmdWithABool("/run/setFastSimulation", this);
    fFastSimulationCmd->SetGuidance("Enable or disable fast simulation.");
    fFastSimulationCmd->SetParameterName("FastSimulation", false);
    fFastSimulationCmd->SetDefaultValue(false);

    fNumberOfScattersCmd = new G4UIcmdWithAnInteger("/run/setNumberOfScatters", this);
    fNumberOfScattersCmd->SetGuidance("Set the number of scatters for fast MC.");
    fNumberOfScattersCmd->SetParameterName("numberOfScattersMax", false);
    fNumberOfScattersCmd->SetDefaultValue(0);

    fMaxEnergyCmd = new G4UIcmdWithADoubleAndUnit("/run/setMaxEnergy", this);
    fMaxEnergyCmd->SetGuidance("Set the maximum energy deposit for fast MC.");    
    fMaxEnergyCmd->SetParameterName("maxEnergy", false);
    fMaxEnergyCmd->SetDefaultValue(0.0);
    fMaxEnergyCmd->SetDefaultUnit("MeV");

    fOutputFileNameCmd = new G4UIcmdWithAString("/run/setOutputFileName", this);
    fOutputFileNameCmd->SetGuidance("Set the output file name.");
    fOutputFileNameCmd->SetParameterName("outputFileName", false);
    fOutputFileNameCmd->SetDefaultValue("G4Sim.root");

}

RunActionMessenger::~RunActionMessenger() {
    delete fFastSimulationCmd;
}

void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
    if (command == fFastSimulationCmd) {
        fRunAction->SetFastSimulation(fFastSimulationCmd->GetNewBoolValue(newValue));
    } else if (command == fNumberOfScattersCmd) {
        fRunAction->SetNumberOfScatters(fNumberOfScattersCmd->GetNewIntValue(newValue));
    } else if (command == fMaxEnergyCmd) {
        fRunAction->SetMaxEnergy(fMaxEnergyCmd->GetNewDoubleValue(newValue));
    } else if (command == fOutputFileNameCmd) {
        fRunAction->SetOutputFileName(newValue);
    }
}

} // namespace G4Sim