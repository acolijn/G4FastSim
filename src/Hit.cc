#include "Hit.hh"
#include "G4SystemOfUnits.hh"


namespace G4FastSim {

G4ThreadLocal G4Allocator<Hit>* HitAllocator = 0;

Hit::Hit()
    : G4VHit(), energyDeposit(0.), position(G4ThreeVector()), time(0.), trackID(-1), parentID(-1), momentum(G4ThreeVector()), particleType(""), processType(""), used(false) {}

Hit::~Hit() {}

void Hit::Print() const {
    G4cout << "TrackID: " << trackID
           << ", ParentID: " << parentID
           << ", Process: " << processType
           << ", Particle: " << particleType
           << ", Position: " << position / cm << " cm"
           << ", Energy Deposit: " << energyDeposit / eV << " eV"
           << ", Time: " << time / ns << " ns"
           << ", Momentum: " << momentum << G4endl;
}

}