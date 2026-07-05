#ifndef SEEBiasingOperator_h
#define SEEBiasingOperator_h 1

#include "G4VBiasingOperator.hh"
#include "G4BOptnChangeCrossSection.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include <map>
#include <atomic>

// Biases the hadronic INELASTIC process for a chosen particle (protons,
// for PANDA) by artificially multiplying its interaction cross section
// by a fixed factor. This follows the standard Geant4 "GB01"-style
// non-splitting cross-section biasing recipe: G4VBiasingOperator
// selects a G4BOptnChangeCrossSection operation for the process to be
// biased; Geant4's biasing framework automatically adjusts the
// interacting track's weight to compensate, so the biased simulation
// remains statistically unbiased in expectation.
//
// IMPORTANT (bug history -- read before touching EventAction/
// SteppingAction scoring): an earlier version of this codebase
// multiplied every accumulated physical quantity (edep, collected
// charge) by the track's weight directly in C++, on the theory that
// "every scored quantity must be weighted". That is WRONG: the charge
// value itself must stay the RAW, unweighted physical quantity for
// that specific simulated history. PANDA_Analyze.py separately
// applies EventWeight when building the histogram (weights=w in
// np.histogram) -- applying weight in BOTH places double-weights every
// biased event, shrinking its charge value by ~1/factor AND its
// histogram count by another ~1/factor, which silently erases the
// entire nuclear-recoil tail once biasing is active. This is invisible
// in a weight=1.0 sanity check (edep*1.0 == edep) and only manifests
// once real fractional weights exist. The correct design: accumulate
// raw edep/charge in EventAction, track the event's representative
// weight separately via EventAction::UpdateEventWeight(), and let
// weight enter the calculation exactly once, at histogram time in
// Python, using the raw charge value for bin placement.
//
// Only the hadronic inelastic process is biased here (see
// ProposeOccurenceBiasingOperation). Elastic scattering and every
// other process for the target particle is left completely unbiased,
// since elastic proton-nucleus recoils in silicon are far too low in
// energy to meaningfully populate the high-charge SEE tail.
class SEEBiasingOperator : public G4VBiasingOperator
{
public:
    // particleToBias: e.g. G4Proton::ProtonDefinition()
    // crossSectionFactor: multiplier applied to the analog (physical)
    //   cross section of the biased process, e.g. 1041 to match
    //   CREME-MC's "Hadronic Cross Section Multiplier". A factor of
    //   1.0 is a no-op and should reproduce the unbiased baseline --
    //   always verify with factor=1.0 before trusting a boosted run.
    SEEBiasingOperator(const G4ParticleDefinition* particleToBias,
                        G4double crossSectionFactor,
                        const G4String& name = "SEEBiasingOperator");
    virtual ~SEEBiasingOperator();

    virtual void StartRun() override;

    // Prints the GLOBAL (all-threads-combined) ground-truth interaction
    // counts to G4cout. Call this explicitly from PANDA.cc right after
    // the run finishes (e.g. after UImanager->ApplyCommand(...) in
    // batch mode), NOT relying on the destructor.
    //
    // BUG HISTORY: an earlier version printed these counts from
    // ~SEEBiasingOperator() and per-instance (non-static) members.
    // That never fired: DetectorConstruction::ConstructSDandField()
    // creates each thread's operator with `new` and never stores or
    // deletes the pointer (a deliberate/common leak pattern in Geant4
    // biasing examples, since operators are meant to live for the
    // run's lifetime) -- so the destructor is simply never called,
    // and the diagnostic silently never printed. Moving to static
    // counters + an explicit call site sidesteps object-lifetime
    // questions entirely.
    static void PrintTotals();

private:
    virtual G4VBiasingOperation*
    ProposeOccurenceBiasingOperation(
        const G4Track* track,
        const G4BiasingProcessInterface* callingProcess) override;

    virtual G4VBiasingOperation*
    ProposeFinalStateBiasingOperation(
        const G4Track*,
        const G4BiasingProcessInterface*) override
    { return nullptr; }

    virtual G4VBiasingOperation*
    ProposeNonPhysicsBiasingOperation(
        const G4Track*,
        const G4BiasingProcessInterface*) override
    { return nullptr; }

    virtual void
    OperationApplied(
        const G4BiasingProcessInterface* callingProcess,
        G4BiasingAppliedCase biasingCase,
        G4VBiasingOperation* occurenceOperationApplied,
        G4double weightForOccurenceInteraction,
        G4VBiasingOperation* finalStateOperationApplied,
        const G4VParticleChange* particleChangeProduced) override;

    const G4ParticleDefinition* fParticleToBias;
    G4double fCrossSectionFactor;

    // Hard counters, independent of any downstream event/CSV scoring.
    // STATIC and ATOMIC: shared across every worker thread's operator
    // instance (Geant4 MT creates one SEEBiasingOperator per thread
    // via ConstructSDandField(), all biasing the same particle/process
    // with the same factor -- these counters give the true combined
    // total across the whole run, read via PrintTotals()).
    //
    // fsNumProposed: how many times ProposeOccurenceBiasingOperation
    //   actually returned a real (non-null) biasing operation for the
    //   target particle's inelastic process, i.e. how many times the
    //   biasing framework was engaged at all for a step.
    // fsNumConfirmedInteractions: how many of those proposed
    //   operations were subsequently confirmed as an ACTUAL
    //   interaction by OperationApplied() (via
    //   SetInteractionOccured()) -- this is the true count of biased
    //   hadronic-inelastic reactions that occurred in this run, ground
    //   truth independent of Recoil_keV/Proton_keV classification
    //   ambiguity in SteppingAction.cc.
    static std::atomic<G4long> fsNumProposed;
    static std::atomic<G4long> fsNumConfirmedInteractions;
    static std::atomic<G4double> fsCrossSectionFactor;

    // One change-cross-section operation per wrapped biasing process
    // interface on the target particle. Populated in StartRun(); at
    // runtime, ProposeOccurenceBiasingOperation() decides per-process
    // whether to actually apply biasing (hadronic inelastic only) or
    // return nullptr to leave that process fully unbiased.
    std::map<const G4BiasingProcessInterface*, G4BOptnChangeCrossSection*>
        fChangeCrossSectionOperations;
};

#endif
