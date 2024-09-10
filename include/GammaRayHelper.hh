#ifndef GAMMA_RAY_HELPER_H
#define GAMMA_RAY_HELPER_H

#include "G4LivermoreComptonModel.hh"
#include "ExtendedLivermoreComptonModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DataVector.hh"
#include "G4AutoLock.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "Hit.hh"
#include <map>
#include <mutex>

/**
 * @namespace G4Sim
 * @brief A namespace for fast simulation classes in Geant4.
/*/
namespace G4Sim {

struct CDFData {
    std::vector<G4double> cdf;
    std::vector<G4double> cosTheta;
};

struct InteractionData {
    G4ThreeVector x0;
    G4ThreeVector dir;
    G4double energy;
    G4double energyDeposited;
    G4double cosTheta;
    G4double weight;
    Hit* hit;
};

/**
 * @class GammaRayHelper
 * @brief A helper class for handling gamma ray interactions.
 *
 * This class provides methods for calculating gamma ray cross sections,
 * Compton scattering angles, and generating interaction points.
 */
class GammaRayHelper {
public:
    static thread_local GammaRayHelper& Instance();
    void Initialize();
    void InitializeCDFs(G4double energy);
    //~GammaRayHelper();


    G4double GetComptonCrossSection(G4double energy, G4Material* material);
    G4double GetPhotoelectricCrossSection(G4double energy, G4Material* material);
    G4double GetTotalCrossSection(G4double energy, G4Material* material);
    G4double GetAttenuationLength(G4double energy, G4Material* material); 
    G4double GetMassAttenuationCoefficient(G4double energy, G4Material* material);

    std::pair<G4ThreeVector, G4double> GenerateInteractionPoint(const G4Step* step);
    std::pair<G4double, G4double> GenerateComptonAngle(const G4Step* step, G4double energy_max);

    InteractionData DoComptonScatter(const G4Step* step, G4ThreeVector x0, G4double energy_max);


    ExtendedLivermoreComptonModel* GetComptonModel(){
        return comptonModel;
    };

    G4LivermorePhotoElectricModel* GetPhotoelectricModel(){
        return photoelectricModel;
    };

private:
    GammaRayHelper();
    ~GammaRayHelper() = default;
    
    CDFData CreateCDF(const G4Element* element, G4double energy);
    G4double CalculateComptonEnergy(G4double energy, G4double cosTheta);
    G4double CalculateMinCosTheta(G4double energy, G4double maxDeposition);
    G4ThreeVector CalculateNewDirection(G4ThreeVector dir, G4double cosTheta);

    GammaRayHelper(const GammaRayHelper&) = delete;
    GammaRayHelper& operator=(const GammaRayHelper&) = delete;

    ExtendedLivermoreComptonModel* comptonModel;
    G4LivermorePhotoElectricModel* photoelectricModel;
    G4LivermoreRayleighModel* rayleighModel;

    std::vector<std::string> fElementsUsed;

    std::map<const G4Element*, CDFData> cdfDataMap; // CDFs for each element

    bool cdfsInitialized;
    G4double fInitialEnergy;


    //std::mutex initMutex; // Mutex for thread-safe initialization
};

}
#endif // GAMMA_RAY_HELPER_H
