#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <atomic>

class G4ParticleGun;
class G4Event;
class G4GenericMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    // /sim/beamXY setter (bound via DeclareMethodWithUnit rather than
    // DeclarePropertyWithUnit) so every worker thread's write also
    // mirrors into the static fsBeamXY below -- same "racy in principle,
    // harmless in practice, every thread writes the same value" pattern
    // SEEBiasingOperator uses for its cross-section factors. Needed so
    // EventAction::PrintUpsetSummary() (a static, thread-independent
    // call site in PANDA.cc) can compute a beam-area cross-section
    // without needing a live PrimaryGeneratorAction instance.
    void SetBeamXY(G4double val);

    static G4double GetBeamXY();

private:
    G4ParticleGun* fParticleGun;
    G4GenericMessenger* fMessenger;

    G4double fEnergy;

    G4double fBeamXY;

    static std::atomic<G4double> fsBeamXY;
};

#endif
