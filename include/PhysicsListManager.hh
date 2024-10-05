#ifndef PHYSICSLISTMANAGER_HH
#define PHYSICSLISTMANAGER_HH

#include "G4VModularPhysicsList.hh"
#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmLivermorePhysics.hh"
//#include "G4OpticalProcessIndex.hh"

class PhysicsListManager {
public:
    PhysicsListManager();  // Constructor
    ~PhysicsListManager(); // Destructor

    G4VModularPhysicsList* CreatePhysicsList(); // Function to set up and return physics list

private:
    bool includeOpticalPhysics=false; // Option to enable optical physics

    // Helper methods to add physics processes
    void AddOpticalPhotonPhysics(G4VModularPhysicsList* physicsList);
};

#endif
