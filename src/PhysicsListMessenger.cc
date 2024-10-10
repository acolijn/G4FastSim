#include "PhysicsListMessenger.hh"
#include "PhysicsListManager.hh"
#include "G4UImanager.hh"

PhysicsListMessenger::PhysicsListMessenger(PhysicsListManager* physicsManager)
    : G4UImessenger(), fPhysicsManager(physicsManager) {
    
    G4UIdirectory* dir = new G4UIdirectory("/physics/");
    dir->SetGuidance("UI commands for the physics setup");
    // Define the command to enable/disable optical physics
    fOpticalCmd = new G4UIcmdWithABool("/physics/setOpticalPhysics", this);
    fOpticalCmd->SetGuidance("Enable or disable optical photon physics");
    fOpticalCmd->SetParameterName("OpticalPhysics", true);
    fOpticalCmd->SetDefaultValue(false);  // Default is false (off)
}

PhysicsListMessenger::~PhysicsListMessenger() {
    delete fOpticalCmd;
}

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
    if (command == fOpticalCmd) {
        // Set the optical physics flag in PhysicsListManager
        G4cout << " PhysicsListMessenger::SetNewValue: Setting optical physics to " << newValue << G4endl;
        fPhysicsManager->SetOpticalPhysics(fOpticalCmd->GetNewBoolValue(newValue));
    }
}
