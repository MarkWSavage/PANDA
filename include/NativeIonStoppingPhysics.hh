#ifndef NativeIonStoppingPhysics_h
#define NativeIonStoppingPhysics_h 1

#include "G4VPhysicsConstructor.hh"

// Overrides G4IonParametrisedLossModel's default ICRU73 stopping-power
// table (registered by G4EmStandardPhysics_option4's ConstructProcess(),
// which runs before this constructor's -- see PANDA.cc registration
// order) to stop silently redirecting any ion with atomic number >= 19
// to iron's (Z=26) scaled table. G4EMLOW8.8 actually ships native
// per-ion ICRU73 data files for every Z from 3 to 80 (confirmed on
// disk: ion_stopping_data/icru73/z<Z>_14.dat), including Au197 (Z=79)
// -- but G4IonDEDXScalingICRU73's hardcoded default range
// (minAtomicNumberIon=19, maxAtomicNumberIon=102) redirects every one
// of those ions to Fe before the file is ever read, regardless of
// whether its own native file exists. Confirmed via
// G4IonDEDXHandler::IsApplicable()/BuildDEDXTable(), which look up
// AtomicNumberBaseIon() -- the (possibly redirected) Z -- before any
// file access.
//
// Fix: register a replacement "ICRU73" table using the same underlying
// G4IonStoppingData source but a G4IonDEDXScalingICRU73 constructed
// with an inverted (empty) range, so AtomicNumberBaseIon() always
// returns the ion's own true Z, never Fe/Ar. G4IonParametrisedLossModel
// ::AddDEDXTable() enforces name uniqueness, so registering ours here
// (before physics tables are built) blocks Geant4's own later default
// registration of the same name silently -- see .cc for the ordering
// argument in full.
//
// For any ion whose native file genuinely doesn't exist (Z>80, or any
// Z/material combination ICRU73 never covered), IsApplicable() simply
// returns false and G4IonParametrisedLossModel falls back to its
// existing internal Bragg/Bethe-Bloch handling exactly as it already
// does today for anything ICRU73 doesn't cover -- this change can only
// ever pick a *better* (native) table for a previously Fe/Ar-redirected
// ion, never remove coverage.
class NativeIonStoppingPhysics : public G4VPhysicsConstructor
{
public:
    NativeIonStoppingPhysics()
    : G4VPhysicsConstructor("NativeIonStoppingPhysics") {}

    void ConstructParticle() override {}
    void ConstructProcess() override;
};

#endif
