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
    fParticleName = "proton";
    fBeamXY = 10*um;

    fMessenger =
        new G4GenericMessenger(this, "/sim/", "Simulation control");

    fMessenger->DeclarePropertyWithUnit("energy", "MeV", fEnergy);
    fMessenger->DeclareProperty("particle", fParticleName);
    fMessenger->DeclarePropertyWithUnit("beamXY", "um", fBeamXY);

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
    auto particle =
        G4ParticleTable::GetParticleTable()->FindParticle(fParticleName);

    fParticleGun->SetParticleDefinition(particle);

    fParticleGun->SetParticleEnergy(fEnergy);

    // Uniform beam over square area
    G4double x =
        (G4UniformRand() - 0.5)*fBeamXY;

    G4double y =
        (G4UniformRand() - 0.5)*fBeamXY;

    // Always start well upstream in vacuum
    // Prevent source from being inside silicon
    auto detector =
        static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction()
        );

    G4double sourceZ =
        -(detector->GetSensitiveThickness()/2.0
           + detector->GetDeadThickness()
           + 1.0*um);

    fParticleGun->SetParticlePosition(
        G4ThreeVector(x,y,sourceZ)
    );

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
