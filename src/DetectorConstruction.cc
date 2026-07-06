#include "DetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4UserLimits.hh"
#include "G4RunManager.hh"
#include "G4ParticleTable.hh"
#include "SEEBiasingOperator.hh"

#include <algorithm>

DetectorConstruction::DetectorConstruction()
{
    fSensitiveLogical = nullptr;

    fSensitiveXY = 10*um;
    fSensitiveThickness = 2*um;
    fDeadThickness = 5*um;
    fStepSize = 0.01*um;

    fMessenger =
        new G4GenericMessenger(this, "/sim/", "Simulation control");

auto& xyCmd =
    fMessenger->DeclarePropertyWithUnit(
        "sensitiveXY",
        "um",
        fSensitiveXY);
xyCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& deadXYCmd =
    fMessenger->DeclarePropertyWithUnit(
        "deadXY",
        "um",
        fDeadXY,
        "Lateral (X/Y) size of the dead layer. Independent of "
        "sensitiveXY -- defaults to the same value (10 um) but does "
        "NOT track it, so set both explicitly if they should match.");
deadXYCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& particleCmd =
    fMessenger->DeclareProperty(
        "particle",
        fParticleName,
        "Primary particle name (e.g. proton, neutron, alpha, deuteron, "
        "triton, He3, GenericIon). Also selects which species' hadronic "
        "inelastic process the SEEBiasingOperator biases -- see "
        "ConstructSDandField() and PANDA.cc.");
particleCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& sensMatCmd =
    fMessenger->DeclareProperty(
        "sensitiveMaterial",
        fSensitiveMaterialName,
        "Sensitive volume material: Si, Ge, GaAs, SiC, or GaN. Also "
        "selects the pair-creation energy and mobility/saturation-"
        "velocity constants SteppingAction uses for the MeV->charge "
        "conversion -- see DetectorConstruction::GetSensitivePairCreation"
        "Energy() and friends for values/sources.");
sensMatCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& deadMatCmd =
    fMessenger->DeclareProperty(
        "deadMaterial",
        fDeadMaterialName,
        "Dead-layer/electrode material: Si, Ge, GaAs, SiC, GaN, SiO2, "
        "Al2O3, or TiO2. Independent of sensitiveMaterial.");
deadMatCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& thickCmd =
    fMessenger->DeclarePropertyWithUnit(
        "sensitiveThickness",
        "um",
        fSensitiveThickness);
thickCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& collectCmd =
    fMessenger->DeclareProperty(
        "useCollectionModel",
        fUseCollectionModel
    );
collectCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& lifetimeCmd =
    fMessenger->DeclarePropertyWithUnit(
        "carrierLifetime",
        "ns",
        fCarrierLifetime
    );
lifetimeCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& fieldCmd =
    fMessenger->DeclarePropertyWithUnit(
        "electricField",
        "volt/m",
        fElectricField
    );

fieldCmd.SetStates(G4State_PreInit, G4State_Idle);


auto& deadCmd =
    fMessenger->DeclarePropertyWithUnit(
        "deadThickness",
        "um",
        fDeadThickness);
deadCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& stepCmd =
    fMessenger->DeclarePropertyWithUnit(
        "stepSize",
        "um",
        fStepSize
    );
stepCmd.SetStates(G4State_PreInit, G4State_Idle);

auto& biasCmd =
    fMessenger->DeclareProperty(
        "biasCrossSectionFactor",
        fBiasCrossSectionFactor,
        "Multiplier on the hadronic inelastic cross section for "
        "whichever particle /sim/particle selects, in the sensitive+"
        "dead volumes. 1.0 = no bias (default). Set to match CREME-MC's "
        "Hadronic Cross Section Multiplier for direct comparability. "
        "ALWAYS verify a run with this left at 1.0 reproduces the "
        "unbiased baseline before trusting a boosted run."
    );
biasCmd.SetStates(G4State_PreInit, G4State_Idle);

}

DetectorConstruction::~DetectorConstruction()
{
    delete fMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto nist = G4NistManager::Instance();

    auto vacuum =
        nist->FindOrBuildMaterial("G4_Galactic");

    auto sensitiveMaterial = ResolveMaterial(fSensitiveMaterialName);
    auto deadMaterial      = ResolveMaterial(fDeadMaterialName);

    // Dynamically size world based on geometry
    G4double worldXY =
        20.0 * std::max(fSensitiveXY, fDeadXY);

    G4double totalThickness =
        fDeadThickness + fSensitiveThickness;

    G4double worldZ =
        50.0 * totalThickness;
    
    if(worldZ < 1*mm)
        worldZ = 1*mm;

    auto solidWorld =
        new G4Box(
            "World",
            worldXY/2,
            worldXY/2,
            worldZ/2
        );

    auto logicWorld =
        new G4LogicalVolume(
            solidWorld,
            vacuum,
            "World"
        );

    auto physWorld =
        new G4PVPlacement(
            nullptr,
            G4ThreeVector(),
            logicWorld,
            "World",
            nullptr,
            false,
            0
        );

    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Dead layer directly below sensitive volume
    auto solidDead =
        new G4Box(
            "DeadLayer",
            fDeadXY/2,
            fDeadXY/2,
            fDeadThickness/2
        );

    auto logicDead =
        new G4LogicalVolume(
            solidDead,
            deadMaterial,
            "DeadLayer"
        );

    fDeadLogical = logicDead;

     auto deadVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // blue
         deadVis->SetForceSolid(true);
         logicDead->SetVisAttributes(deadVis);

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(
            0,
            0,
            -(fSensitiveThickness/2 + fDeadThickness/2)
        ),
        logicDead,
        "DeadLayer",
        logicWorld,
        false,
        0
    );

    // Sensitive volume
    auto solidSensitive =
        new G4Box(
            "SensitiveVolume",
            fSensitiveXY/2,
            fSensitiveXY/2,
            fSensitiveThickness/2
        );


    fSensitiveLogical =
        new G4LogicalVolume(
            solidSensitive,
            sensitiveMaterial,
            "SensitiveVolume"
        );

    auto sensVis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
        sensVis->SetForceSolid(true);
        fSensitiveLogical->SetVisAttributes(sensVis);

    fSensitiveLogical->SetUserLimits(
        //new G4UserLimits(0.01*um)     //Appropriate for nano-electronics
        //new G4UserLimits(0.50*um)     //Appropriate for large structures
        new G4UserLimits(fStepSize)
    );

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0,0,0),
        fSensitiveLogical,
        "SensitiveVolume",
        logicWorld,
        false,
        0
    );



    return physWorld;
}

G4LogicalVolume* DetectorConstruction::GetSensitiveLogical() const
{
    return fSensitiveLogical;
}

G4Material* DetectorConstruction::ResolveMaterial(const G4String& name)
{
    auto nist = G4NistManager::Instance();

    if (name == "Si")    return nist->FindOrBuildMaterial("G4_Si");
    if (name == "Ge")    return nist->FindOrBuildMaterial("G4_Ge");
    if (name == "GaAs")  return nist->FindOrBuildMaterial("G4_GALLIUM_ARSENIDE");
    if (name == "SiO2")  return nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    if (name == "Al2O3") return nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
    if (name == "TiO2")  return nist->FindOrBuildMaterial("G4_TITANIUM_DIOXIDE");

    // SiC and GaN aren't in Geant4's NIST compound database -- build
    // them from elements. G4NistManager::FindOrBuildMaterial caches by
    // name internally, but "SiC"/"GaN" aren't NIST names it recognizes,
    // so we cache these ourselves to avoid rebuilding (and re-warning
    // Geant4's material table about a duplicate name) on every call.
    if (name == "SiC")
    {
        static G4Material* siliconCarbide = nullptr;
        if (!siliconCarbide)
        {
            siliconCarbide = new G4Material("SiC", 3.21*g/cm3, 2);
            siliconCarbide->AddElement(nist->FindOrBuildElement("Si"), 1);
            siliconCarbide->AddElement(nist->FindOrBuildElement("C"), 1);
        }
        return siliconCarbide;
    }

    if (name == "GaN")
    {
        static G4Material* galliumNitride = nullptr;
        if (!galliumNitride)
        {
            galliumNitride = new G4Material("GaN", 6.15*g/cm3, 2);
            galliumNitride->AddElement(nist->FindOrBuildElement("Ga"), 1);
            galliumNitride->AddElement(nist->FindOrBuildElement("N"), 1);
        }
        return galliumNitride;
    }

    G4Exception(
        "DetectorConstruction::ResolveMaterial()",
        "InvalidMaterial",
        FatalException,
        ("Unknown material name: '" + name +
         "' -- expected one of Si, Ge, GaAs, SiC, GaN, SiO2, Al2O3, "
         "TiO2.").c_str()
    );
    return nullptr;
}

G4double DetectorConstruction::GetSensitivePairCreationEnergy() const
{
    // Mean energy to create one electron-hole pair. Typical room-
    // temperature bulk literature values -- treat as approximate,
    // same spirit as the "(tune these)" mobility/vsat comment below.
    if (fSensitiveMaterialName == "Si")   return 3.6*eV;
    if (fSensitiveMaterialName == "Ge")   return 2.9*eV;
    if (fSensitiveMaterialName == "GaAs") return 4.2*eV;
    if (fSensitiveMaterialName == "SiC")  return 7.6*eV;
    if (fSensitiveMaterialName == "GaN")  return 8.9*eV;

    G4Exception(
        "DetectorConstruction::GetSensitivePairCreationEnergy()",
        "InvalidMaterial",
        FatalException,
        ("Unknown /sim/sensitiveMaterial name: '" + fSensitiveMaterialName +
         "' -- expected one of Si, Ge, GaAs, SiC, GaN.").c_str()
    );
    return 3.6*eV;
}

G4double DetectorConstruction::GetSensitiveElectronMobility() const
{
    // Low-field electron mobility. Typical room-temperature bulk
    // literature values (tune these).
    if (fSensitiveMaterialName == "Si")   return 1350.0 * cm2/volt/second;
    if (fSensitiveMaterialName == "Ge")   return 3900.0 * cm2/volt/second;
    if (fSensitiveMaterialName == "GaAs") return 8500.0 * cm2/volt/second;
    if (fSensitiveMaterialName == "SiC")  return  900.0 * cm2/volt/second;
    if (fSensitiveMaterialName == "GaN")  return 1200.0 * cm2/volt/second;

    G4Exception(
        "DetectorConstruction::GetSensitiveElectronMobility()",
        "InvalidMaterial",
        FatalException,
        ("Unknown /sim/sensitiveMaterial name: '" + fSensitiveMaterialName +
         "' -- expected one of Si, Ge, GaAs, SiC, GaN.").c_str()
    );
    return 1350.0 * cm2/volt/second;
}

G4double DetectorConstruction::GetSensitiveSaturationVelocity() const
{
    // Typical room-temperature bulk literature values (tune these).
    // GaAs's electron velocity-field curve is non-monotonic (Gunn
    // effect / negative differential mobility): it peaks around
    // 2.0e7 cm/s near ~3-4 kV/cm, then droops before re-settling at
    // high field. This simple mu*E-capped-at-vsat model can't
    // reproduce that NDR dip, so 2.0e7 (the commonly tabulated
    // "saturation velocity" figure, e.g. Sze & Ng) is used as the cap
    // -- more realistic than reusing Si's 1.0e7, which was an
    // oversimplification from the first cut of this model.
    if (fSensitiveMaterialName == "Si")   return 1.0e7 * cm/second;
    if (fSensitiveMaterialName == "Ge")   return 6.0e6 * cm/second;
    if (fSensitiveMaterialName == "GaAs") return 2.0e7 * cm/second;
    if (fSensitiveMaterialName == "SiC")  return 2.0e7 * cm/second;
    if (fSensitiveMaterialName == "GaN")  return 2.5e7 * cm/second;

    G4Exception(
        "DetectorConstruction::GetSensitiveSaturationVelocity()",
        "InvalidMaterial",
        FatalException,
        ("Unknown /sim/sensitiveMaterial name: '" + fSensitiveMaterialName +
         "' -- expected one of Si, Ge, GaAs, SiC, GaN.").c_str()
    );
    return 1.0e7 * cm/second;
}

void DetectorConstruction::ConstructSDandField()
{
    // Attach the SEE biasing operator to the sensitive and dead-layer
    // volumes. This is called once (serial) or once per worker thread
    // (MT), after Construct() -- the Geant4-recommended, thread-safe
    // place to attach biasing operators (each thread gets its own
    // operator instance, matching Geant4's official biasing examples).
    //
    // With fBiasCrossSectionFactor == 1.0 (the default), the operator
    // still runs but proposes no actual change to the physical cross
    // section, so results should be statistically identical to not
    // having biasing installed at all. This is the required sanity
    // check before trusting any boosted-factor run -- see the
    // messenger command help text and SEEBiasingOperator.hh for why.
    //
    // The particle to bias is whatever /sim/particle selected (read
    // via fParticleName, set before /run/initialize triggers this
    // method). PANDA.cc must have wrapped that same particle's
    // hadronic inelastic process via PhysicsBias() at physics-
    // construction time -- if it hasn't (e.g. a particle name not in
    // that list), SEEBiasingOperator::StartRun() prints an explicit
    // WARNING and biasing is silently inert for this run.
    G4ParticleDefinition* primaryDef =
        G4ParticleTable::GetParticleTable()->FindParticle(fParticleName);

    if (!primaryDef)
    {
        G4Exception(
            "DetectorConstruction::ConstructSDandField()",
            "InvalidParticle",
            FatalException,
            ("Unknown /sim/particle name: '" + fParticleName +
             "' -- not found in G4ParticleTable.").c_str()
        );
    }

    auto* biasingOperator =
        new SEEBiasingOperator(
            primaryDef,
            fBiasCrossSectionFactor
        );

    if (fSensitiveLogical)
        biasingOperator->AttachTo(fSensitiveLogical);

    if (fDeadLogical)
        biasingOperator->AttachTo(fDeadLogical);
}
