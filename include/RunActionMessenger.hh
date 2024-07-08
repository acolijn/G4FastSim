#ifndef RUN_ACTION_MESSENGER_H
#define RUN_ACTION_MESSENGER_H 1

#include "G4UImessenger.hh"
#include "RunAction.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"


namespace G4FastSim {

class RunAction;

class RunActionMessenger : public G4UImessenger {
public:
    RunActionMessenger(RunAction* action);
    ~RunActionMessenger() override;

    void SetNewValue(G4UIcommand* command, G4String newValue) override;

private:
    RunAction* fRunAction;

    G4UIcmdWithABool* fFastSimulationCmd;
    G4UIcmdWithAnInteger* fNumberOfScattersCmd;

};

} // namespace G4FastSim

#endif
