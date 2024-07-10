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

//GammaRayHelper::GammaRayHelper() {}
GammaRayHelper::GammaRayHelper() : comptonModel(nullptr), photoelectricModel(nullptr), cdfsInitialized(false) {}


/**
 * Initialize the Compton and photoelectric models.
 */
void GammaRayHelper::Initialize() {

    G4cout << "GammaRayHelper::Initializing Compton and photoelectric models " << G4endl;

    comptonModel = new ExtendedLivermoreComptonModel();
    photoelectricModel = new G4LivermorePhotoElectricModel();

    G4DataVector cuts;
    cuts.push_back(1 * keV); // Example cut value, adjust as needed?? what it means? find out.....
    comptonModel->Initialise(G4Gamma::Gamma(), cuts);
    photoelectricModel->Initialise(G4Gamma::Gamma(), cuts);
} 

/**
 * @brief Initializes the cumulative distribution functions (CDFs) for gamma rays.
 * 
 * This function initializes the CDFs for gamma rays based on the given energy.
 * It checks if the CDFs have already been initialized and if not, performs the initialization.
 * The CDFs are created for each element used in the material table.
 * 
 * @param energy The energy of the gamma rays.
 */
void GammaRayHelper::InitializeCDFs(G4double energy) {

    fInitialEnergy = energy;

    if (!cdfsInitialized) {
        Initialize();
        fElementsUsed.clear();
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        for (size_t i = 0; i < materialTable->size(); ++i) {
            G4Material* material = (*materialTable)[i];
            const G4ElementVector* elementVector = material->GetElementVector();
            for (size_t j = 0; j < material->GetNumberOfElements(); ++j) {
                const G4Element* elm = (*elementVector)[j];
                if (std::find(fElementsUsed.begin(), fElementsUsed.end(), elm->GetName()) == fElementsUsed.end()) {
                    fElementsUsed.push_back(elm->GetName());
                    cdfDataMap[elm] = CreateCDF(elm, energy); // Create initial CDF for a typical energy
                }
            }
        }
        cdfsInitialized = true;
    }
}

/**
 * Create a cumulative distribution function (CDF) for a given element and energy.
 * @param element The element to create the CDF for.
 * @param energy The energy of the gamma ray.
 * @return The CDF data.
 */
CDFData GammaRayHelper::CreateCDF(const G4Element* element, G4double energy) {
    G4int nPoints = 1000;
    G4double cosThetaMin = -1.0;
    G4double cosThetaMax =  1.0;
    G4double dCosTheta = (cosThetaMax - cosThetaMin) / nPoints;

    std::vector<G4double> cdf(nPoints + 1);
    std::vector<G4double> cosTheta(nPoints + 1);

    G4double sum = 0.0;
    for (G4int i = 0; i <= nPoints; ++i) {
        G4double cosThetaVal = cosThetaMin + i * dCosTheta;
        cosTheta[i] = cosThetaVal;
        G4double diffCrossSection = comptonModel->DifferentialCrossSection(element, new G4DynamicParticle(G4Gamma::Gamma(), G4ThreeVector(1,0,0), energy), cosThetaVal) * dCosTheta;
        sum += diffCrossSection;
        cdf[i] = sum;
    }

    for (G4int i = 0; i <= nPoints; ++i) {
        cdf[i] /= sum;
    }

    CDFData cdfData;
    cdfData.cdf = cdf;
    cdfData.cosTheta = cosTheta;

    return cdfData;
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

G4double GammaRayHelper::GenerateComptonAngle(const G4Step* step) {

    const G4MaterialCutsCouple* couple = step->GetPreStepPoint()->GetMaterialCutsCouple();
    G4double energy0 = step->GetPreStepPoint()->GetKineticEnergy();
    const G4ParticleDefinition* particle = G4Gamma::Gamma();
    const G4Element* element = comptonModel->SelectRandomAtom(couple, particle, energy0);

    // Check if a new CDF needs to be generated
    CDFData cdfData;
    if (energy0 != fInitialEnergy) {
        cdfData = CreateCDF(element, energy0);
    } else {
        cdfData = cdfDataMap[element];
    }

    G4double rand = G4UniformRand();
    const std::vector<G4double>& cdf = cdfData.cdf;
    const std::vector<G4double>& cosTheta = cdfData.cosTheta;

    auto it = std::lower_bound(cdf.begin(), cdf.end(), rand);
    G4int idx = std::distance(cdf.begin(), it);

    G4double cosThetaVal = cosTheta.back();
    if (idx > 0 && idx < cdf.size()) {
        G4double t1 = cdf[idx - 1];
        G4double t2 = cdf[idx];
        G4double cosTheta1 = cosTheta[idx - 1];
        G4double cosTheta2 = cosTheta[idx];
        cosThetaVal = cosTheta1 + (cosTheta2 - cosTheta1) * (rand - t1) / (t2 - t1);
    }

    return cosThetaVal;
}

}
 