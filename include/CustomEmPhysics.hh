#ifndef CustomEmPhysics_hh
#define CustomEmPhysics_hh 1

#include "G4VPhysicsConstructor.hh"
#include "G4EmLivermorePhysics.hh"

class CustomEmPhysics : public G4EmLivermorePhysics {
public:
    CustomEmPhysics();
    virtual ~CustomEmPhysics();

    virtual void ConstructProcess();
};

#endif
