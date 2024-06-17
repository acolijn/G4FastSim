#include "GammaRayHelper.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4AutoLock.hh"
#include <thread>

// Declare the mutex globally within the file for G4AutoLock
namespace {
    G4Mutex mutex = G4MUTEX_INITIALIZER;
}

GammaRayHelper& GammaRayHelper::Instance() {
    static thread_local GammaRayHelper instance;
    return instance;
}

GammaRayHelper::GammaRayHelper() {}

/* GammaRayHelper::~GammaRayHelper() {
    // Clean up allocated memory
    for (auto& pair : comptonModels) {
        delete pair.second;
    }
    for (auto& pair : photoelectricModels) {
        delete pair.second;
    }
} */

void GammaRayHelper::Initialize(G4Material* mat) {
    G4AutoLock lock(&mutex); // Use G4AutoLock for thread-safe initialization
    if (comptonModels.find(mat) == comptonModels.end()) {
        G4cout << "GammaRayHelper::Initializing Compton and photoelectric models for material " << mat->GetName() << G4endl;
        G4LivermoreComptonModel* comptonModel = new G4LivermoreComptonModel();
        G4LivermorePhotoElectricModel* photoelectricModel = new G4LivermorePhotoElectricModel();

        G4DataVector cuts;
        cuts.push_back(1.0 * keV); // Example cut value, adjust as needed
        comptonModel->Initialise(G4Gamma::Gamma(), cuts);
        photoelectricModel->Initialise(G4Gamma::Gamma(), cuts);

        comptonModels[mat] = comptonModel;
        photoelectricModels[mat] = photoelectricModel;
    }
} 
 
G4ThreeVector GammaRayHelper::GenerateComptonScatteringDirection(
    const G4ThreeVector& initialDirection,
    G4double initialEnergy,
    G4double& scatteredEnergy,
    G4double minAngle,
    G4double maxAngle,
    G4double& weight,
    G4Material* material)
{
    auto& comptonModel = comptonModels[material];

    // Calculate the scattering angle and corresponding differential cross section
    G4double theta, phi;
    G4double minCosTheta = std::cos(maxAngle);
    G4double maxCosTheta = std::cos(minAngle);

    // Ensure correct min/max range
    if (minCosTheta > maxCosTheta) std::swap(minCosTheta, maxCosTheta);

    G4double cosTheta, sinTheta, differentialCrossSection;
    do {
        cosTheta = minCosTheta + (maxCosTheta - minCosTheta) * G4UniformRand();
        sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
        phi = 2.0 * pi * G4UniformRand();
        
        differentialCrossSection = 0.0;
        const G4ElementVector* elementVector = material->GetElementVector();
        const G4double* fractionVector = material->GetFractionVector();
        for (size_t i = 0; i < material->GetNumberOfElements(); ++i) {
            const G4Element* element = (*elementVector)[i];
            differentialCrossSection += fractionVector[i] * comptonModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), initialEnergy, element->GetZ(), cosTheta);
        }
    } while (G4UniformRand() > differentialCrossSection);

    theta = std::acos(cosTheta);

    // Calculate the energy of the scattered photon using Compton formula
    G4double scatteredWavelength = (h_Planck * c_light / initialEnergy) + (h_Planck / (electron_mass_c2 * c_light)) * (1.0 - cosTheta);
    scatteredEnergy = h_Planck * c_light / scatteredWavelength;

    // Calculate the new direction of the scattered photon
    G4double sinPhi = std::sin(phi);
    G4double cosPhi = std::cos(phi);

    G4ThreeVector newDirection(
        sinTheta * cosPhi,
        sinTheta * sinPhi,
        cosTheta
    );

    // Rotate the new direction to align with the initial photon direction
    G4ThreeVector rotationAxis = G4ThreeVector(0., 0., 1.).cross(initialDirection.unit());
    G4double rotationAngle = std::acos(G4ThreeVector(0., 0., 1.).dot(initialDirection.unit()));
    G4RotationMatrix rotationMatrix;
    rotationMatrix.rotate(rotationAngle, rotationAxis);

    newDirection *= rotationMatrix;

    // Calculate the weight
    G4double fullCrossSection = 0.0;
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* fractionVector = material->GetFractionVector();
    for (size_t i = 0; i < material->GetNumberOfElements(); ++i) {
        const G4Element* element = (*elementVector)[i];
        fullCrossSection += fractionVector[i] * comptonModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), initialEnergy, element->GetZ(), -1.0); // Total cross section
    }
    G4double restrictedCrossSection = (maxCosTheta - minCosTheta) * fullCrossSection / 2.0;
    weight = fullCrossSection / restrictedCrossSection;

    return newDirection;
}

G4double GammaRayHelper::GetComptonCrossSection(G4double energy, G4Material* material) {
    //G4AutoLock lock(&mutex); // Ensure thread-safe access

    
    auto& comptonModel = comptonModels[material];
    G4double crossSection = 0.0;
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* fractionVector = material->GetFractionVector();
    for (size_t i = 0; i < material->GetNumberOfElements(); ++i) {
        const G4Element* element = (*elementVector)[i];
        crossSection += fractionVector[i] * comptonModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ());
    }
    return crossSection;
}

G4double GammaRayHelper::GetPhotoelectricCrossSection(G4double energy, G4Material* material) {
    //G4AutoLock lock(&mutex); // Ensure thread-safe access

    auto& photoelectricModel = photoelectricModels[material];
    G4double crossSection = 0.0;
    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* fractionVector = material->GetFractionVector();
    for (size_t i = 0; i < material->GetNumberOfElements(); ++i) {
        const G4Element* element = (*elementVector)[i];
        crossSection += fractionVector[i] * photoelectricModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ());
    }
    return crossSection;
}


G4double GammaRayHelper::GetTotalCrossSection(G4double energy, G4Material* material) {
    //G4AutoLock lock(&mutex); // Ensure thread-safe access
    G4double crossSection = GetComptonCrossSection(energy, material) + GetPhotoelectricCrossSection(energy, material);
    return crossSection;
}

G4double GammaRayHelper::GetAttenuationLength(G4double energy, G4Material* material) {
    //G4AutoLock lock(&mutex); // Ensure thread-safe access
    
    G4double crossSection = GetComptonCrossSection(energy, material) + GetPhotoelectricCrossSection(energy, material);
    return 1.0 / (crossSection * material->GetDensity());
}

G4double GammaRayHelper::GetMassAttenuationCoefficient(G4double energy, G4Material* material) {
    /*
    
    This function calculates the mass attenuation coefficient for a given material at a given energy.
 

    A.P. Colijn, 2024
    */
    G4double att = 0; 

    auto& photoelectricModel = photoelectricModels[material];
    auto& comptonModel = comptonModels[material];

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
 