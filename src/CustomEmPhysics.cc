#include "CustomEmPhysics.hh"
#include "G4ProcessManager.hh"
#include "G4RayleighScattering.hh"
#include "G4Gamma.hh"

CustomEmPhysics::CustomEmPhysics()
    : G4EmLivermorePhysics() {}

CustomEmPhysics::~CustomEmPhysics() {}

void CustomEmPhysics::ConstructProcess() {
    // Call the base class method to construct the standard processes
    G4EmLivermorePhysics::ConstructProcess();

    // Get the process manager for gamma particles
    G4ProcessManager* pManager = G4Gamma::Gamma()->GetProcessManager();

    // Find and remove Rayleigh scattering
    G4RayleighScattering* rayleighScattering = nullptr;
    G4int nProcesses = pManager->GetProcessListLength();
    for (G4int i = 0; i < nProcesses; ++i) {
        G4VProcess* process = (*pManager->GetProcessList())[i];
        if ((rayleighScattering = dynamic_cast<G4RayleighScattering*>(process)) != nullptr) {
            G4cout << "CustomEmPhysics::ConstructProcess:  Removing Rayleigh scattering" << G4endl;	
            pManager->RemoveProcess(rayleighScattering);
            break;
        }
    }
}
