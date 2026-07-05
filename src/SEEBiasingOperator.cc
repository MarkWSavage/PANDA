#include "SEEBiasingOperator.hh"

#include "G4BiasingProcessInterface.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessType.hh"
#include "G4HadronicProcessType.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "G4Threading.hh"

#include <cfloat>

// Static member definitions (shared across every worker thread's
// operator instance -- see SEEBiasingOperator.hh for why).
std::atomic<G4long> SEEBiasingOperator::fsNumProposed{0};
std::atomic<G4long> SEEBiasingOperator::fsNumConfirmedInteractions{0};
std::atomic<G4double> SEEBiasingOperator::fsCrossSectionFactor{1.0};

SEEBiasingOperator::SEEBiasingOperator(
    const G4ParticleDefinition* particleToBias,
    G4double crossSectionFactor,
    const G4String& name)
: G4VBiasingOperator(name),
  fParticleToBias(particleToBias),
  fCrossSectionFactor(crossSectionFactor)
{
    // All threads construct their operator with the same factor (it
    // comes from the same run.mac setting), so this write is racy in
    // principle but harmless in practice -- every thread writes the
    // same value.
    fsCrossSectionFactor.store(crossSectionFactor);
}

SEEBiasingOperator::~SEEBiasingOperator()
{
    // NOTE: do not rely on this destructor for diagnostic output --
    // see PrintTotals() and the bug-history comment in the header for
    // why. This destructor may never run in practice (the operator is
    // created with `new` in DetectorConstruction::ConstructSDandField()
    // and never explicitly deleted, matching common Geant4 biasing
    // example patterns).
    for (auto& kv : fChangeCrossSectionOperations)
        delete kv.second;
}

void SEEBiasingOperator::PrintTotals()
{
    G4long proposed = fsNumProposed.load();
    G4long confirmed = fsNumConfirmedInteractions.load();
    G4double factor = fsCrossSectionFactor.load();

    G4cout << G4endl;
    G4cout << "=======================================================" << G4endl;
    G4cout << "SEEBiasingOperator -- GLOBAL FINAL COUNTS" << G4endl;
    G4cout << "(all worker threads combined; ground truth, independent"
           << " of CSV/Recoil_keV/Proton_keV scoring)" << G4endl;
    G4cout << "    Cross section factor used   : " << factor << G4endl;
    G4cout << "    Proposed biasing operations : " << proposed << G4endl;
    G4cout << "    Confirmed real interactions : " << confirmed << G4endl;
    if (proposed > 0)
    {
        G4cout << "    Confirmed/Proposed ratio    : "
               << (G4double)confirmed / proposed
               << " (expect ~= 1/crossSectionFactor = "
               << 1.0 / factor << ")" << G4endl;
    }
    else
    {
        G4cout << "    WARNING: zero proposed operations -- biasing "
               << "never engaged at all. Check StartRun() output above "
               << "for 'wrapped process(es)' confirmation." << G4endl;
    }
    G4cout << "=======================================================" << G4endl;
    G4cout << G4endl;
}

void SEEBiasingOperator::StartRun()
{
    // Create one change-cross-section operation for every physics
    // process wrapped for biasing on the target particle (this is
    // populated by G4GenericBiasingPhysics::PhysicsBias(...) called
    // in PANDA.cc, which explicitly wraps only "protonInelastic").
    // Each is a candidate; ProposeOccurenceBiasingOperation below
    // decides at runtime whether a given wrapped process is actually
    // the hadronic inelastic process we want to bias -- all others
    // are left unbiased (nullptr returned for them).
    const G4ProcessManager* processManager =
        fParticleToBias->GetProcessManager();

    const G4BiasingProcessSharedData* sharedData =
        G4BiasingProcessInterface::GetSharedData(processManager);

    if (!sharedData)
    {
        G4cout << "SEEBiasingOperator::StartRun() -- WARNING: no biasing "
               << "shared data found for " << fParticleToBias->GetParticleName()
               << ". Biasing will have NO EFFECT. Check that PhysicsBias(...) "
               << "was called for this particle in PANDA.cc." << G4endl;
        return;
    }

    G4cout << "SEEBiasingOperator::StartRun() -- wrapped process(es) for "
           << fParticleToBias->GetParticleName() << " biasing:" << G4endl;

    for (std::size_t i = 0;
         i < sharedData->GetPhysicsBiasingProcessInterfaces().size();
         ++i)
    {
        const G4BiasingProcessInterface* wrapperProcess =
            sharedData->GetPhysicsBiasingProcessInterfaces()[i];

        G4String wrappedName =
            wrapperProcess->GetWrappedProcess()->GetProcessName();

        G4cout << "    - " << wrappedName << G4endl;

        G4String operationName = "XSboost-" + wrappedName;

        fChangeCrossSectionOperations[wrapperProcess] =
            new G4BOptnChangeCrossSection(operationName);
    }
}

G4VBiasingOperation* SEEBiasingOperator::ProposeOccurenceBiasingOperation(
    const G4Track* track,
    const G4BiasingProcessInterface* callingProcess)
{
    if (track->GetDefinition() != fParticleToBias)
        return nullptr;

    // Only bias the hadronic INELASTIC process. Elastic scattering
    // and every other process for this particle is left completely
    // unbiased -- returning nullptr restores normal, analog physics
    // for them.
    const G4VProcess* wrapped = callingProcess->GetWrappedProcess();

    if (wrapped->GetProcessType() != fHadronic)
        return nullptr;

    if (wrapped->GetProcessSubType() != fHadronInelastic)
        return nullptr;

    G4double analogInteractionLength =
        wrapped->GetCurrentInteractionLength();

    // Process not active for this step/track (e.g. below threshold,
    // or particle about to leave the biased volume): leave it alone.
    if (analogInteractionLength > DBL_MAX / 10.0)
        return nullptr;

    G4double analogXS = 1.0 / analogInteractionLength;

    G4BOptnChangeCrossSection* operation =
        fChangeCrossSectionOperations[callingProcess];

    if (!operation)
        return nullptr;

    operation->SetBiasedCrossSection(fCrossSectionFactor * analogXS);
    operation->Sample();

    ++fsNumProposed;

    return operation;
}

void SEEBiasingOperator::OperationApplied(
    const G4BiasingProcessInterface* callingProcess,
    G4BiasingAppliedCase /*biasingCase*/,
    G4VBiasingOperation* occurenceOperationApplied,
    G4double /*weightForOccurenceInteraction*/,
    G4VBiasingOperation* /*finalStateOperationApplied*/,
    const G4VParticleChange* /*particleChangeProduced*/)
{
    G4BOptnChangeCrossSection* operation =
        fChangeCrossSectionOperations[callingProcess];

    if (operation == occurenceOperationApplied)
    {
        operation->SetInteractionOccured();
        ++fsNumConfirmedInteractions;
    }
}
