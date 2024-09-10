#include "Materials.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4SystemOfUnits.hh"


namespace G4Sim {

Materials::Materials()
{
}

Materials::~Materials()
{
}

void Materials::DefineMaterials()
{
    G4cout << "Materials::DefineMaterials" << G4endl;

    // Get nist material manager
    G4NistManager* nist = G4NistManager::Instance();
    //==== Elements ====
    G4cout<<"Materials::DefineMaterials: Elements"<<G4endl;
    G4Element *H = nist->FindOrBuildElement("H");
    G4Element *C = nist->FindOrBuildElement("C");
    G4Element *N = nist->FindOrBuildElement("N");
    G4Element *O = nist->FindOrBuildElement("O");
    G4Element *F = nist->FindOrBuildElement("F");
    G4Element *Si = nist->FindOrBuildElement("Si");
    G4Element *Fe = nist->FindOrBuildElement("Fe");
    G4Element *Cr = nist->FindOrBuildElement("Cr");
    G4Element *Ni = nist->FindOrBuildElement("Ni");
    G4Element *Xe = nist->FindOrBuildElement("Xe");


    //==== Materials ====

    //==== Air ====
    G4cout << "Materials::DefineMaterials: Air" << G4endl;
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
    //==== Water ====
    G4cout << "Materials::DefineMaterials: Water" << G4endl;
    G4Material* water = nist->FindOrBuildMaterial("G4_WATER");
    //==== Liquid Xenon ====
    G4Material *LXe = new G4Material("LXe", 2.862 * g / cm3, 1, kStateLiquid,
                                       177.05 * kelvin, 1.5 * atmosphere);
    // DR 20180518 - Density according to:
    // -
    // https://xe1t-wiki.lngs.infn.it/doku.php?id=xenon:xenon1t:deg:tpc:targetmass
    // -
    // https://xe1t-wiki.lngs.infn.it/doku.php?id=xenon:xenon1t:analysis:sciencerun1:sc_summary
    LXe->AddElement(Xe, 1);

    //==== Vacuum ====
    G4cout << "Materials::DefineMaterials: Vacuum" << G4endl;
    G4Material *Vacuum = new G4Material("Vacuum", 1.e-25 * g / cm3, 2, kStateGas);
    Vacuum->AddElement(N, 0.755);
    Vacuum->AddElement(O, 0.245);
    //==== Teflon / PTFE ====
    G4cout << "Materials::DefineMaterials: PTFE" << G4endl;
    G4Material *PTFE = new G4Material("PTFE", 2.2 * g / cm3, 2, kStateSolid);
    PTFE->AddElement(C, 0.240183);
    PTFE->AddElement(F, 0.759817);
    //==== Stainless Steel ====
    G4cout << "Materials::DefineMaterials: StainlessSteel" << G4endl;
    G4Material *StainlessSteel = new G4Material("StainlessSteel", 8.06 * g / cm3, 3, kStateSolid);
    StainlessSteel->AddElement(Fe, 0.70);
    StainlessSteel->AddElement(Cr, 0.18);
    StainlessSteel->AddElement(Ni, 0.12);

    G4cout << "Materials::DefineMaterials done" << G4endl;
}

} // namespace G4Sim