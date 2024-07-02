#include "RunActionMessenger.hh"
#include "RunAction.hh"
#include "G4UIcmdWithABool.hh"

///using namespace G4FastSim;

RunActionMessenger::RunActionMessenger(RunAction* action)
    : fRunAction(action) {

    fFastSimulationCmd = new G4UIcmdWithABool("/run/setFastSimulation", this);
    fFastSimulationCmd->SetGuidance("Enable or disable fast simulation.");
    fFastSimulationCmd->SetParameterName("FastSimulation", false);
    fFastSimulationCmd->SetDefaultValue(false);
}

RunActionMessenger::~RunActionMessenger() {
    delete fFastSimulationCmd;
}

void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
    if (command == fFastSimulationCmd) {
        fRunAction->SetFastSimulation(fFastSimulationCmd->GetNewBoolValue(newValue));
    }
}

//} // namespace G4FastSim