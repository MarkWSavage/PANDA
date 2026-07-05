#include "DetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4UserLimits.hh"
#include "G4RunManager.hh"
#include "G4Proton.hh"
#include "SEEBiasingOperator.hh"

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
        "Multiplier on the hadronic inelastic cross section for protons "
        "in the sensitive+dead volumes. 1.0 = no bias (default). Set to "
        "match CREME-MC's Hadronic Cross Section Multiplier for direct "
        "comparability. ALWAYS verify a run with this left at 1.0 "
        "reproduces the unbiased baseline before trusting a boosted run."
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

    auto silicon =
        nist->FindOrBuildMaterial("G4_Si");

    // Dynamically size world based on geometry
    G4double worldXY = 20.0 * fSensitiveXY;

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
            fSensitiveXY/2,
            fSensitiveXY/2,
            fDeadThickness/2
        );

    auto logicDead =
        new G4LogicalVolume(
            solidDead,
            silicon,
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
            silicon,
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
    auto* biasingOperator =
        new SEEBiasingOperator(
            G4Proton::ProtonDefinition(),
            fBiasCrossSectionFactor
        );

    if (fSensitiveLogical)
        biasingOperator->AttachTo(fSensitiveLogical);

    if (fDeadLogical)
        biasingOperator->AttachTo(fDeadLogical);
}
