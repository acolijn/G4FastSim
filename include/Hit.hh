#ifndef HIT_HH
#define HIT_HH

#include "G4VHit.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
#include "G4String.hh"
#include "G4THitsCollection.hh"

namespace G4FastSim {

class Hit : public G4VHit {
public:
    Hit();
    virtual ~Hit();

    // Hit properties
    G4double energyDeposit;
    G4ThreeVector position;
    G4double time;
    G4int trackID;
    G4int parentID;
    G4ThreeVector momentum;
    G4String particleType;

    // Operators
    inline void* operator new(size_t);
    inline void operator delete(void* hit);

private:
};

extern G4ThreadLocal G4Allocator<Hit>* HitAllocator;

inline void* Hit::operator new(size_t) {
    if (!HitAllocator) HitAllocator = new G4Allocator<Hit>;
    return (void*)HitAllocator->MallocSingle();
}

inline void Hit::operator delete(void* hit) {
    HitAllocator->FreeSingle((Hit*)hit);
}

// Define the HitsCollection type
typedef G4THitsCollection<Hit> HitsCollection;


}

#endif
