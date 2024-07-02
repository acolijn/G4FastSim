#include "Hit.hh"

namespace G4FastSim {

G4ThreadLocal G4Allocator<Hit>* HitAllocator = 0;

Hit::Hit()
    : G4VHit(), energyDeposit(0.), position(G4ThreeVector()), time(0.), trackID(-1), parentID(-1), momentum(G4ThreeVector()), particleType("") {}

Hit::~Hit() {}

}
