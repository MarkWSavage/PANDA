#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fParticleGun = new G4ParticleGun(1);

    // Defaults
    fEnergy = 50*MeV;
    fBeamXY = 10*um;

    fMessenger =
        new G4GenericMessenger(this, "/sim/", "Simulation control");

    fMessenger->DeclarePropertyWithUnit("energy", "MeV", fEnergy);
    fMessenger->DeclarePropertyWithUnit("beamXY", "um", fBeamXY);
    // /sim/particle is declared by DetectorConstruction, not here --
    // see DetectorConstruction::GetParticleName() for why.

    // Beam direction toward +Z
    fParticleGun->SetParticleMomentumDirection(
        G4ThreeVector(0,0,1)
    );
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
    delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // Particle name lives on DetectorConstruction (see
    // DetectorConstruction::GetParticleName() for why) since the
    // biasing operator also needs it.
    auto detector =
        static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction()
        );

    auto particle =
        G4ParticleTable::GetParticleTable()
            ->FindParticle(detector->GetParticleName());

    fParticleGun->SetParticleDefinition(particle);

    fParticleGun->SetParticleEnergy(fEnergy);

    // Uniform beam over square area
    G4double x =
        (G4UniformRand() - 0.5)*fBeamXY;

    G4double y =
        (G4UniformRand() - 0.5)*fBeamXY;

    // Always start well upstream in vacuum
    // Prevent source from being inside silicon (or, if a lid is
    // present, inside/behind the lid -- see DetectorConstruction's
    // fLidThickness. GetLidThickness()/GetLidGap() are 0 when no lid
    // is configured, so this collapses back to the original offset.)

    G4double sourceZ =
        -(detector->GetSensitiveThickness()/2.0
           + detector->GetDeadThickness()
           + detector->GetLidGap()
           + detector->GetLidThickness()
           + 1.0*um);

    fParticleGun->SetParticlePosition(
        G4ThreeVector(x,y,sourceZ)
    );

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
