#ifndef DetectorConstruction_h
#define DetectorConstruction_h
#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4GenericMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    // Constructs and attaches the SEE biasing operator to the
    // sensitive and dead-layer logical volumes. Called once (serial)
    // or once per worker thread (MT) by the run manager, after
    // Construct(). This is the Geant4-recommended, thread-safe place
    // to attach biasing operators -- see SEEBiasingOperator.hh.
    virtual void ConstructSDandField();

    G4LogicalVolume* GetSensitiveLogical() const;
    G4LogicalVolume* GetDeadLogical() const { return fDeadLogical; }
    G4double GetSensitiveThickness() const { return fSensitiveThickness; }
    G4double GetDeadThickness() const { return fDeadThickness; }
    G4double GetSensitiveXY() const { return fSensitiveXY; }
    G4bool GetUseCollectionModel() const { return fUseCollectionModel; }
    G4double GetCarrierLifetime() const { return fCarrierLifetime; }
    G4double GetElectricField() const { return fElectricField; }
    G4double GetBiasCrossSectionFactor() const { return fBiasCrossSectionFactor; }
private:
    G4LogicalVolume* fSensitiveLogical;
    G4LogicalVolume* fDeadLogical = nullptr;
    G4GenericMessenger* fMessenger;

    G4double fStepSize = 0.018*um;
    G4double fSensitiveXY = 10*um;
    G4double fSensitiveThickness = 10*um;
    G4double fDeadThickness = 5*um;

    G4bool   fUseCollectionModel = true;
    G4double fCarrierLifetime    = 10.0 * ns;
    G4double fElectricField      = 1.0 * kilovolt/cm;

    // Multiplier applied to the hadronic inelastic cross section for
    // protons in the sensitive + dead volumes (see SEEBiasingOperator).
    // 1.0 = no bias (default); set via /sim/biasCrossSectionFactor.
    // ALWAYS verify a factor=1.0 run reproduces the unbiased baseline
    // before trusting results from a boosted factor.
    G4double fBiasCrossSectionFactor = 1.0;

};

#endif
