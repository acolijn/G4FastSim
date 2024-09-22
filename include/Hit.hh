#ifndef HIT_HH
#define HIT_HH

#include "G4VHit.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
#include "G4String.hh"
#include "G4THitsCollection.hh"

/**
 * @namespace G4Sim
 * @brief A namespace for fast simulation classes in Geant4.
/*/
namespace G4Sim {

/**
 * @brief The Hit class represents a hit in the detector.
 * 
 * It inherits from the G4VHit class and contains properties such as energy deposit, position, time, track ID, parent ID,
 * momentum, particle type, and process type.
 */
class Hit : public G4VHit {
public:
    /**
     * @brief Default constructor for the Hit class.
     */
    Hit();

    /**
     * @brief Destructor for the Hit class.
     */
    virtual ~Hit();

    // Hit properties
    G4double energyDeposit; /**< Energy deposit of the hit. */
    G4ThreeVector position; /**< Position of the hit. */
    G4double time; /**< Time of the hit. */
    G4int trackID; /**< Track ID of the hit. */
    G4int parentID; /**< Parent ID of the hit. */
    G4ThreeVector momentum; /**< Momentum of the hit. */
    G4String particleType; /**< Particle type of the hit. */
    G4String processType; /**< Process type of the hit. */
    G4double particleEnergy0; /**< Energy of the particle at the beginning of a step */
    G4double particleEnergy1; /**< Energy of the particle after the step */
    G4bool used; /**< Flag to indicate if the hit has been used in a cluster. */

    // New: Tag for identifying hit collection (default: logical volume name)
    G4String tag;
    
    /**
     * @brief Print function to print the details of the hit.
     */
    void Print() const;

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
