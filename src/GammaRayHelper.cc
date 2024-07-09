#include "GammaRayHelper.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4AutoLock.hh"
#include "G4DynamicParticle.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4LivermoreComptonModel.hh"
#include <thread>
#include <vector>

// Declare the mutex globally within the file for G4AutoLock
namespace {
    G4Mutex mutex = G4MUTEX_INITIALIZER;
}

namespace G4FastSim {

GammaRayHelper& GammaRayHelper::Instance() {
    static thread_local GammaRayHelper instance;
    return instance;
}

GammaRayHelper::GammaRayHelper() {}

/**
 * Initialize the Compton and photoelectric models.
 */
void GammaRayHelper::Initialize() {

    G4cout << "GammaRayHelper::Initializing Compton and photoelectric models " << G4endl;

    comptonModel = new ExtendedLivermoreComptonModel();
    photoelectricModel = new G4LivermorePhotoElectricModel();

    G4DataVector cuts;
    cuts.push_back(1.0 * keV); // Example cut value, adjust as needed
    comptonModel->Initialise(G4Gamma::Gamma(), cuts);
    photoelectricModel->Initialise(G4Gamma::Gamma(), cuts);

} 

/**
 * Get the Compton cross section for a given energy and material.
 * @param energy The energy of the gamma ray.
 * @param material The material to calculate the cross section for.
 * @return The Compton cross section.
 */
G4double GammaRayHelper::GetComptonCrossSection(G4double energy, G4Material* material) {
    //G4AutoLock lock(&mutex); // Ensure thread-safe access
    G4double crossSection = 0.0;
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* fractionVector = material->GetFractionVector();
    for (size_t i = 0; i < material->GetNumberOfElements(); ++i) {
        const G4Element* element = (*elementVector)[i];
        crossSection += fractionVector[i] * comptonModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ());
    }
    return crossSection;
}

/**
 * Get the photoelectric cross section for a given energy and material.
 * @param energy The energy of the gamma ray.
 * @param material The material to calculate the cross section for.
 * @return The photoelectric cross section.
 */
G4double GammaRayHelper::GetPhotoelectricCrossSection(G4double energy, G4Material* material) {
    //G4AutoLock lock(&mutex); // Ensure thread-safe access
    G4double crossSection = 0.0;
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* fractionVector = material->GetFractionVector();
    for (size_t i = 0; i < material->GetNumberOfElements(); ++i) {
        const G4Element* element = (*elementVector)[i];
        crossSection += fractionVector[i] * photoelectricModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ());
    }
    return crossSection;
}

/**
 * Get the total cross section (Compton + photoelectric) for a given energy and material.
 * @param energy The energy of the gamma ray.
 * @param material The material to calculate the cross section for.
 * @return The total cross section.
 */
G4double GammaRayHelper::GetTotalCrossSection(G4double energy, G4Material* material) {
    //G4AutoLock lock(&mutex); // Ensure thread-safe access
    G4double crossSection = GetComptonCrossSection(energy, material) + GetPhotoelectricCrossSection(energy, material);
    return crossSection;
}

/**
 * Get the attenuation length for a given energy and material.
 * @param energy The energy of the gamma ray.
 * @param material The material to calculate the attenuation length for.
 * @return The attenuation length.
 */
G4double GammaRayHelper::GetAttenuationLength(G4double energy, G4Material* material) {
    //G4AutoLock lock(&mutex); // Ensure thread-safe access
    
    G4double crossSection = GetComptonCrossSection(energy, material) + GetPhotoelectricCrossSection(energy, material);
    return 1.0 / (crossSection * material->GetDensity());
}

/**
 * Get the mass attenuation coefficient for a given energy and material.
 * @param energy The energy of the gamma ray.
 * @param material The material to calculate the mass attenuation coefficient for.
 * @return The mass attenuation coefficient.
 */
G4double GammaRayHelper::GetMassAttenuationCoefficient(G4double energy, G4Material* material) {
    /*
    
    This function calculates the mass attenuation coefficient for a given material at a given energy.
 

    A.P. Colijn, 2024
    */
    G4double att = 0; 

    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* fractionVector = material->GetFractionVector();

    for (size_t i = 0; i < material->GetNumberOfElements(); ++i) {
        const G4Element* element = (*elementVector)[i];
        G4double A = element->GetAtomicMassAmu();

        G4double sigma = photoelectricModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ());
        sigma += comptonModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ());
        // not 100% sure about this conversion...... maybe I should check teh units of teh cross section routines.....
        att += fractionVector[i] * (sigma / A ) * cm2;
    }
    
    return att;
}

/**
 * Generate the direction of a Compton-scattered gamma ray.
 * @param initialDirection The initial direction of the gamma ray.
 * @param initialEnergy The initial energy of the gamma ray.
 * @param scatteredEnergy The energy of the scattered gamma ray (output parameter).
 * @param minAngle The minimum scattering angle.
 * @param maxAngle The maximum scattering angle.
 * @param weight The weight of the scattered gamma ray (output parameter).
 * @param material The material in which the scattering occurs.
 * @param step The step in which the scattering occurs.
 * @return The direction of the scattered gamma ray.
 */
G4ThreeVector GammaRayHelper::GenerateComptonScatteringDirection(
    //const G4ThreeVector& initialDirection,
    //G4double initialEnergy,
    //G4double& scatteredEnergy,
    //G4double minAngle,
    //G4double maxAngle,
    //G4double& weight,
    G4Material* material, const G4Step *step)
{

    // first select the material to which you want to scatter......   

    G4cout << "GammaRayHelper::GenerateComptonScatteringDirection" << G4endl; 

    const G4MaterialCutsCouple* couple = step->GetPreStepPoint()->GetMaterialCutsCouple();
    G4double energy0 = 1.1 * MeV;
    const G4ParticleDefinition* particle = G4Gamma::Gamma();
    const G4Element *element = GetComptonModel()->SelectRandomAtom(couple, particle, energy0);

    //   G4cout << "Selected element: " << element->GetName() << G4endl;

    //G4DynamicParticle* gamma = new G4DynamicParticle(G4Gamma::Gamma(), initialDirection, initialEnergy);
    //comptonModel->SampleSecondaries(fv, cuts, gamma, 0,0);

    //G4DynamicParticle* electron = (*fv)[0];
    //scatteredEnergy = electron->GetKineticEnergy();

    G4ThreeVector newDirection = G4ThreeVector(0,0,1.);
    return newDirection;
}

}
 