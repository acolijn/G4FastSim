#ifndef PHYSICSLISTMESSENGER_HH
#define PHYSICSLISTMESSENGER_HH

#include "G4UImessenger.hh"
#include "G4UIcmdWithABool.hh"

class PhysicsListManager;

class PhysicsListMessenger : public G4UImessenger {
public:
    PhysicsListMessenger(PhysicsListManager* physicsManager);
    virtual ~PhysicsListMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValue) override;

private:
    PhysicsListManager* fPhysicsManager;  // Pointer to PhysicsListManager

    G4UIcmdWithABool* fOpticalCmd;  // Command to enable/disable optical physics
};

#endif