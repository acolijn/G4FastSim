#include "PhysicsListManager.hh"
#include "PhysicsListMessenger.hh"

PhysicsListManager::PhysicsListManager() { 
    // make a messenger
    fMessenger = new PhysicsListMessenger(this);
}

PhysicsListManager::~PhysicsListManager() { }

G4VModularPhysicsList* PhysicsListManager::CreatePhysicsList() {
    G4PhysListFactory factory;
    
    // Choose a reference physics list as the base
    G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("FTFP_BERT_HP");
    // Replace the default EM physics with Livermore for better precision in low-energy EM interactions
    physicsList->ReplacePhysics(new G4EmLivermorePhysics());

    // Add optical photon physics if enabled
    if (includeOpticalPhysics) {
        G4cout << "PhysicsListManager::CreatePhysicsList: Adding optical photon physics" << G4endl;
        AddOpticalPhotonPhysics(physicsList);
    } else {
        G4cout << "PhysicsListManager::CreatePhysicsList: Optical photon physics is disabled" << G4endl;
    }

    return physicsList;
}

void PhysicsListManager::AddOpticalPhotonPhysics(G4VModularPhysicsList* physicsList) {
    // Create an optical physics module
    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
    physicsList->RegisterPhysics(opticalPhysics);
}
