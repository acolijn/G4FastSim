#ifndef GAMMA_RAY_HELPER_H
#define GAMMA_RAY_HELPER_H

#include "G4LivermoreComptonModel.hh"
#include "ExtendedLivermoreComptonModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DataVector.hh"
#include "G4AutoLock.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include <map>
#include <mutex>

namespace G4FastSim {


class GammaRayHelper {
public:
    static thread_local GammaRayHelper& Instance();
    void Initialize();
    //~GammaRayHelper();

    G4ThreeVector GenerateComptonScatteringDirection(
        //const G4ThreeVector& initialDirection,
        //G4double initialEnergy,
        //G4double& scatteredEnergy,
        //G4double minAngle,
        //G4double maxAngle,
        //G4double& weight,
        G4Material* material, const G4Step *step);

    G4double GetComptonCrossSection(G4double energy, G4Material* material);
    G4double GetPhotoelectricCrossSection(G4double energy, G4Material* material);
    G4double GetTotalCrossSection(G4double energy, G4Material* material);
    G4double GetAttenuationLength(G4double energy, G4Material* material); 
    G4double GetMassAttenuationCoefficient(G4double energy, G4Material* material);

    ExtendedLivermoreComptonModel* GetComptonModel(){
        return comptonModel;
    };

    G4LivermorePhotoElectricModel* GetPhotoelectricModel(){
        return photoelectricModel;
    };

private:
    GammaRayHelper();
    ~GammaRayHelper() = default;

    GammaRayHelper(const GammaRayHelper&) = delete;
    GammaRayHelper& operator=(const GammaRayHelper&) = delete;

    ExtendedLivermoreComptonModel* comptonModel;
    G4LivermorePhotoElectricModel* photoelectricModel;
    //std::mutex initMutex; // Mutex for thread-safe initialization
};

}
#endif // GAMMA_RAY_HELPER_H
