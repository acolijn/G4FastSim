#include "Hit.hh"
#include "G4SystemOfUnits.hh"


namespace G4Sim {

G4ThreadLocal G4Allocator<Hit>* HitAllocator = 0;

Hit::Hit()
    : G4VHit(), energyDeposit(0.), position(G4ThreeVector()), time(0.), trackID(-1), parentID(-1), momentum(G4ThreeVector()), particleType(""), processType(""), used(false), tag("") {}

Hit::~Hit() {}

void Hit::Print() const {
    G4cout << "ID: " << trackID
           << ", ID: " << parentID
           << ", Proc: " << processType
           << ", Ptcl: " << particleType
           << ", E0: " << particleEnergy0 / keV << " keV"
           << ", E1: " << particleEnergy1 / keV << " keV"
           << ", dE: " << energyDeposit / keV << " keV"
           << ", pos: " << position / cm << " cm"
           << ", t: " << time / ns << " ns"
           << G4endl;
    //       << ", p: " << momentum << G4endl;
}

}