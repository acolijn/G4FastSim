#include "ExtendedLivermoreComptonModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4Gamma.hh"
#include "G4LogLogInterpolation.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

/**
 * @brief Constructs an instance of the ExtendedLivermoreComptonModel class.
 *
 * This constructor initializes the ExtendedLivermoreComptonModel object by calling the base class constructor G4LivermoreComptonModel().
 * It also creates a scatterInterpolation object of type G4LogLogInterpolation and sets the scatterFile to "comp/ce-sf-".
 * The scatterFunctionData is then initialized as a new G4CompositeEMDataSet object with scatterInterpolation and load the data from scatterFile.
 */
ExtendedLivermoreComptonModel::ExtendedLivermoreComptonModel() : G4LivermoreComptonModel() {
    G4VDataSetAlgorithm* scatterInterpolation = new G4LogLogInterpolation;
    G4String scatterFile = "comp/ce-sf-";
    scatterFunctionData = new G4CompositeEMDataSet(scatterInterpolation, 1., 1.);
    scatterFunctionData->LoadData(scatterFile);
}

ExtendedLivermoreComptonModel::~ExtendedLivermoreComptonModel() {

    delete scatterFunctionData;
}

/**
 * Calculates the scatter function for Compton scattering in the Extended Livermore model.
 * 
 * @param elm The element to calculate the scatter function for.
 * @param aDynamicGamma The dynamic gamma particle.
 * @param theta The scattering angle.
 * @return The scatter function value.
 */
G4double ExtendedLivermoreComptonModel::FormFactor(const G4Element* elm, const G4DynamicParticle* aDynamicGamma, G4double theta) {
    // Get the energy of the gamma
    G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();

    // Get the particle 
    const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
    // Get Z of the element
    G4int Z = (G4int)elm->GetZ();

    // Calculate the wavelength of the photon
    G4double wlPhoton = h_Planck*c_light/photonEnergy0;
    G4double oneCost = 1-std::cos(theta);
    G4double x = std::sqrt(oneCost/2.) / (wlPhoton/cm);

    // return the scatter function value
    G4double f = scatterFunctionData->FindValue(x,Z-1);

    return f;
}

/**
 * Calculates the Klein-Nishina scattering function value for a given gamma particle and scattering angle.
 *
 * @param aDynamicGamma The dynamic gamma particle.
 * @param theta The scattering angle in radians.
 * @return The Klein-Nishina scattering function value.
 */
G4double ExtendedLivermoreComptonModel::KleinNishina(const G4DynamicParticle* aDynamicGamma, G4double theta) {

    // Get the energy of the gamma
    G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();
    // Calculate the wavelength of the photon
    G4double cost = std::cos(theta);
    G4double e0m = photonEnergy0/electron_mass_c2;

    G4double ff = 1 + e0m*(1-cost);
    G4double kn = (1. + cost*cost + e0m*e0m*(1.-cost)*(1.-cost)/ff)/(ff*ff)/2.;

    // return the scatter function value
    return kn;
}