#include "GammaRayHelper.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4AutoLock.hh"
#include "G4DynamicParticle.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include <thread>
#include <vector>

// Declare the mutex globally within the file for G4AutoLock
namespace {
    G4Mutex mutex = G4MUTEX_INITIALIZER;
}

namespace G4Sim {

GammaRayHelper& GammaRayHelper::Instance() {
    static thread_local GammaRayHelper instance;
    return instance;
}

//GammaRayHelper::GammaRayHelper() {}
GammaRayHelper::GammaRayHelper() : comptonModel(nullptr), photoelectricModel(nullptr), rayleighModel(nullptr), cdfsInitialized(false) {}


/**
 * Initialize the Compton and photoelectric models.
 */
void GammaRayHelper::Initialize() {

    G4cout << "GammaRayHelper::Initializing Compton and photoelectric models " << G4endl;

    comptonModel = new ExtendedLivermoreComptonModel();
    photoelectricModel = new G4LivermorePhotoElectricModel();
    rayleighModel = new G4LivermoreRayleighModel();

    G4DataVector cuts;
    cuts.push_back(0 * keV); // Example cut value, adjust as needed?? what it means? find out.....

    comptonModel->Initialise(G4Gamma::Gamma(), cuts);
    photoelectricModel->Initialise(G4Gamma::Gamma(), cuts);
    rayleighModel->Initialise(G4Gamma::Gamma(), cuts);
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
    G4int nPoints = 10000;
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

    const G4ElementVector* elementVector = material->GetElementVector();
    const G4double* fractionVector = material->GetFractionVector();
    G4double mass_attenuation = 0.0;
    for (size_t i = 0; i < material->GetNumberOfElements(); ++i) {
        const G4Element* element = (*elementVector)[i];
        //G4double xsec = photoelectricModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ()) + 
        //                comptonModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ()) +
        //                rayleighModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ());
        G4double xsec = photoelectricModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ()) + 
                        comptonModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ());
        G4double A = element->GetA();
        //G4cout << i << " Element: " << element->GetName() << " xsec: " << xsec << " A: " << A << G4endl;
        //G4cout << i << "             phot  = " << comptonModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ())  << G4endl;
        //G4cout << i << "             compt = " << comptonModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ())  << G4endl;
        //G4cout << i << "             rayl  = " << rayleighModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ())  << G4endl;
        //G4cout << i<<" Element: " << element->GetName() << " xsec: " << xsec << " A: " << A << " r ="<< comptonModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(), energy, element->GetZ())/xsec<<G4endl;
        mass_attenuation += fractionVector[i] * xsec * Avogadro / A ;
        
    }
    G4double rho = material->GetDensity();
    G4double attenuationLength = 1/(rho*mass_attenuation);

    //G4cout << "Material: " << material->GetName() << G4endl;
    //G4cout << "Energy: " << energy/keV << " keV"<<G4endl;
    //G4cout << "Density: " << material->GetDensity() / (g / cm3)<< " g/cm3"<< G4endl;
    //G4cout << "Attenuation length: " << attenuationLength / cm<< " cm"<<G4endl;

    return attenuationLength;
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
 * Generates a Compton scattering angle for a gamma ray.
 * @param step The G4Step object containing information about the current step.
 * @param energy_max The maximum energy deposition allowed for the gamma ray.
 * @return The cosine of the scattering angle.
 */
std::pair<G4double,G4double> GammaRayHelper::GenerateComptonAngle(const G4Step* step, G4double energy_max) {

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

    // get the minimal value of cos(theta)
    G4double cosThetaMin = CalculateMinCosTheta(energy0, energy_max);
    // get corresponding CDF value..
    G4double cdfMin = 0;
    auto it0 = std::lower_bound(cdfData.cosTheta.begin(), cdfData.cosTheta.end(), cosThetaMin);
    if (it0 != cdfData.cosTheta.end()) {
        G4int idx = std::distance(cdfData.cosTheta.begin(), it0);
        cdfMin = cdfData.cdf[idx];
    }
    G4double weight = 1 - cdfMin;

    // Generate a random number between cdfMin and 1
    G4double rand = G4UniformRand() * (1 - cdfMin) + cdfMin;
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

  return std::make_pair(cosThetaVal,weight);
}

/**
 * Calculates the energy of the scattered gamma ray after a Compton scattering event.
 *
 * @param energy0 The initial energy of the gamma ray.
 * @param cosTheta The scattering angle of the gamma ray.
 * @return The energy of the scattered gamma ray.
 */
G4double GammaRayHelper::CalculateComptonEnergy(G4double energy0, G4double cosTheta) {
    G4double energy1 = energy0 / (1 + energy0 / (511 * keV) * (1 - cosTheta));
    return energy1; 
}

/**
 * Calculates the minimum cosine of the scattering angle for a gamma ray.
 *
 * @param energy0 The initial energy of the gamma ray.
 * @param maxDeposition The maximum energy deposition allowed for the gamma ray.
 * @return The minimum cosine of the scattering angle.
 */
G4double GammaRayHelper::CalculateMinCosTheta(G4double energy0, G4double maxDeposition) {
    // if maxDeposition is larger than the energy, the minimum cosTheta is -1
    if (maxDeposition >= energy0) {
        return -1;
    }
    G4double minCosTheta = 1.0 - 511 * keV * maxDeposition / energy0 / (energy0 - maxDeposition);

    if (minCosTheta > 1) {
        minCosTheta = 1;
    } else if (minCosTheta < -1) {
        minCosTheta = -1;
    }

    return minCosTheta;
}

/**
 * Calculates the new direction of the particle after a Compton scatter.
 *
 * @param step The G4Step object containing information about the current step.
 * @param cosTheta The cosine of the scattering angle.
 * @return The new direction of the particle.
 */
G4ThreeVector GammaRayHelper::CalculateNewDirection(G4ThreeVector direction, G4double cosTheta) {
  
  G4double phi = twopi * G4UniformRand();
  G4double cosPhi = std::cos(phi);
  G4double sinPhi = std::sin(phi);
  G4double sinTheta = std::sqrt(1 - cosTheta * cosTheta);

  G4ThreeVector newDirection(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
  // now the magic happens. The new direction rotated towards the original direction
  G4ThreeVector rotationAxis = G4ThreeVector(0, 0, 1).cross(direction);
  G4double rotationAngle = std::acos(G4ThreeVector(0, 0, 1).dot(direction));

  G4RotationMatrix rotationMatrix;
  rotationMatrix.rotate(rotationAngle, rotationAxis);
  // Transform the new direction back to the original coordinate system
  newDirection = rotationMatrix * newDirection;

  return newDirection;
}

InteractionData GammaRayHelper::DoComptonScatter(const G4Step* step, G4ThreeVector x0, G4double energy_max) {
    // Get the gamma ray energy
    G4double energy0 = step->GetPreStepPoint()->GetKineticEnergy();
    // Generate the Compton scattering angle
    auto result = GenerateComptonAngle(step, energy_max);
    G4double cosTheta = result.first;
    G4double weight = result.second;
    // Calculate the energy of the scattered gamma ray
    G4double energy1 = CalculateComptonEnergy(energy0, cosTheta);
    // Calculate the new direction of the gamma ray
    G4ThreeVector direction = CalculateNewDirection(step->GetPreStepPoint()->GetMomentumDirection(), cosTheta);
    // 
    // // G4cout << "Compton scattering: " << energy0/keV << " -> " << energy1/keV << " cos(theta) = " << cosTheta << " weight = " << weight << " emax = "<<energy_max/keV<<G4endl;

    InteractionData data;
    data.x0 = x0;
    data.dir = direction;
    data.energy = energy1;
    data.energyDeposited = energy0 - energy1;
    data.cosTheta = cosTheta;
    data.weight = weight;
    data.hit = new Hit();

    return data;
}


/**
 * Generates an interaction point along a line segment between two points.
 *
 * @param entrance The entrance point of the line segment.
 * @param exit The exit point of the line segment.
 * @param attenuation_length The attenuation length of the material.
 * @return The generated interaction point.
 */
std::pair<G4ThreeVector, G4double> GammaRayHelper::GenerateInteractionPoint(const G4Step *step) {

  // get the energy of the particle
  G4double energy = step->GetPreStepPoint()->GetKineticEnergy();
  
  G4LogicalVolume* volume_pre = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4Material* material        = volume_pre->GetMaterial();
  G4double attenuation_length = GetAttenuationLength(energy, material);

  // maximum possible distance....
  G4double maxDistance = step->GetStepLength();

  // get the entrance and exit points of the step
  G4ThreeVector entrance = step->GetPreStepPoint()->GetPosition();
  G4ThreeVector exit = step->GetPostStepPoint()->GetPosition();

  // generate a random distance along the line segment
  G4double rand = G4UniformRand();
  // We want to generate a hit inside teh fiducial volume. 
  // Therefore there is a maximum to the CDF that we can generate: this is returned for weighing the event
  G4double maxCDF = 1 - std::exp(-maxDistance / attenuation_length);
  G4double adjustedRand = rand * maxCDF;
  G4double distance = -attenuation_length * std::log(1 - adjustedRand);
  G4ThreeVector interactionPoint = entrance + (exit - entrance).unit() * distance;

  return std::make_pair(interactionPoint,maxCDF);
}



}