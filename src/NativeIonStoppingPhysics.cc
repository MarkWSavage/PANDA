#include "NativeIonStoppingPhysics.hh"

#include "G4GenericIon.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4IonStoppingData.hh"
#include "G4IonDEDXScalingICRU73.hh"
#include "G4EmParameters.hh"
#include "G4SystemOfUnits.hh"

void NativeIonStoppingPhysics::ConstructProcess()
{
    // GenericIon's "ionIoni" process/model is built by
    // G4EmStandardPhysics_option4::ConstructProcess() (registered in
    // PANDA.cc before this constructor via ReplacePhysics), which
    // G4VModularPhysicsList guarantees runs its ConstructProcess()
    // before ours -- registration order is call order for both
    // ConstructParticle() and ConstructProcess() across all physics
    // constructors on every thread (master and each MT worker alike),
    // so this is safe under PANDA's multithreaded run manager with no
    // separate per-worker hook needed.
    G4ParticleDefinition* ion = G4GenericIon::GenericIon();
    G4ProcessManager* pManager = ion->GetProcessManager();

    if (!pManager)
        return;

    G4ProcessVector* pv = pManager->GetProcessList();

    for (G4int i = 0; i < pv->size(); ++i)
    {
        G4VProcess* proc = (*pv)[i];

        if (proc->GetProcessName() != "ionIoni")
            continue;

        auto* eLossProc = dynamic_cast<G4VEnergyLossProcess*>(proc);
        if (!eLossProc)
            continue;

        // Only one model is registered for GenericIon's ionIoni process
        // when G4EmStandardPhysics_option4 is constructed with a
        // non-empty DEDX name (PANDA.cc passes "ICRU73") -- see
        // G4EmStandardPhysics_option4::ConstructProcess(), which calls
        // ionIoni->SetEmModel(new G4IonParametrisedLossModel()) at
        // index 0 in that case, with no separate low/high-energy model
        // split for ions (unlike protons/alphas).
        auto* model =
            dynamic_cast<G4IonParametrisedLossModel*>(eLossProc->EmModel(0));

        if (!model)
            continue;

        G4bool icru90 = G4EmParameters::Instance()->UseICRU90Data();

        // Same table name ("ICRU73") that
        // G4IonParametrisedLossModel::Initialise() uses for its own
        // default registration later, during physics-table building
        // (triggered by /run/initialize, well after every physics
        // constructor's ConstructProcess() has already run). AddDEDXTable
        // enforces name uniqueness and silently refuses a second table
        // under a name already present -- registering ours first here
        // means Geant4's own later default attempt is the one that's
        // rejected, not ours.
        // Scoped, not blanket: validated against SRIM across PANDA's
        // heavy-ion roster (C12/Ni58/I127/Au197), switching to native
        // per-ion data for EVERY Z>=19 made Ni58 (Z=28, only 2 above
        // Fe's 26) worse -- the existing Fe-scaling approximation is
        // already good for an ion this close to the reference, and the
        // native file for it is not obviously better validated. I127
        // (Z=53, 27 away) came out roughly a wash (flipped from +6.4%
        // over to -5.6% under). Au197 (Z=79, 53 away) improved clearly
        // (+17.8% -> +8.4%). Cutoff set at the midpoint between I127 and
        // Au197's Z (66): Fe-scaling still applies for 19<=Z<=66
        // (unchanged behavior for Ni58/I127 and anything in that band),
        // native data only for Z>66 (Au197 and, if ever used, any
        // similarly heavy or heavier primary up to Z=80 where native
        // files exist -- beyond that, IsApplicable() finding no cached
        // table for the missing Z falls back to the model's existing
        // Bragg/Bethe-Bloch handling exactly as any ICRU73-uncovered
        // case already does today).
        model->AddDEDXTable(
            "ICRU73",
            new G4IonStoppingData("ion_stopping_data/icru", icru90),
            new G4IonDEDXScalingICRU73(19, 66));

        break;
    }
}
