#ifndef __EXTENDED_LIVERMORE_COMPTON_MODEL_H__
#define __EXTENDED_LIVERMORE_COMPTON_MODEL_H__

#include "G4LivermoreComptonModel.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DynamicParticle.hh"
#include "G4VEMDataSet.hh"

class ExtendedLivermoreComptonModel : public G4LivermoreComptonModel {
    public:
        ExtendedLivermoreComptonModel();
        virtual ~ExtendedLivermoreComptonModel();

        G4double FormFactor(const G4Element* , const G4DynamicParticle* aDynamicGamma, G4double cosTheta);
        G4double KleinNishina(const G4DynamicParticle* aDynamicGamma, G4double cosTheta);
        G4double DifferentialCrossSection(const G4Element* element, const G4DynamicParticle* aDynamicGamma, G4double cossTheta);

    private:
        G4VEMDataSet* scatterFunctionData;
};

#endif // __EXTENDED_LIVERMORE_COMPTON_MODEL_H__