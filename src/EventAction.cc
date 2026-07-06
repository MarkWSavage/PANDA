#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "EventAction.hh"
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4ios.hh"

namespace {
    const char* kResultsDir = "Results/Current";

    // Geant4 MT's master thread has G4GetThreadId() == -1; a plain
    // (non-MT) G4RunManager never spawns worker threads at all, so
    // this constructor only ever runs with threadID == -1 in that
    // case -- either way, "events_t-1.csv" is a valid, unique
    // filename, so no special-casing is needed for serial vs MT.
    std::string ThreadEventsPath()
    {
        return std::string(kResultsDir) + "/events_t"
             + std::to_string(G4Threading::G4GetThreadId()) + ".csv";
    }
}

EventAction::EventAction()
: G4UserEventAction()
{
    fMessenger = new G4GenericMessenger(this, "/sim/", "Simulation control");

    fMessenger->DeclareMethod(
        "criticalCharge",
        &EventAction::SetCriticalCharge,
        "Set critical charge in fC");

    fMessenger->DeclareProperty(
        "verbose",
        fVerbose,
        "Enable verbose output");

    std::filesystem::create_directories(kResultsDir);

    fCSV.open(ThreadEventsPath());

    fCSV << "EventID,"
         << "DepositedCharge_fC,"
         << "CollectedCharge_fC,"
         << "Total_keV,"
         << "Proton_keV,"
         << "Electron_keV,"
         << "Recoil_keV,"
         << "UpsetCharge_fC,"
         << "EventWeight"
         << G4endl;
}

EventAction::~EventAction()
{
    if (fCSV.is_open())
        fCSV.close();

    delete fMessenger;
}

void EventAction::BeginOfEventAction(const G4Event*)
{
    fHits.clear();

    fTotalEdep    = 0.0;
    fProtonEdep   = 0.0;
    fElectronEdep = 0.0;
    fRecoilEdep   = 0.0;
    fCollectedCharge = 0.0;

    fEventWeight = 1.0;
    fMaxSingleEdep = 0.0;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    // Charge equivalent of the RAW deposited energy (ideal case: every
    // e-h pair fully collected, no trapping/recombination loss).
    G4double depositedCharge =
        (fTotalEdep / (3.6 * eV)) * CLHEP::eplus;

    G4double depositedCharge_fC =
        depositedCharge / (1.0e-15 * CLHEP::coulomb);

    // Charge equivalent ACTUALLY collected at the electrode, after the
    // drift/trapping model in SteppingAction (always computed there,
    // regardless of which one is used below for the upset criterion).
    G4double collectedCharge_fC =
        fCollectedCharge / (1.0e-15 * CLHEP::coulomb);

    auto detector =
        static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()
                ->GetUserDetectorConstruction()
        );

    // Which of the two charges is used to decide whether this event
    // causes an upset is controlled by /sim/useCollectionModel.
    G4double upsetCharge =
        detector->GetUseCollectionModel()
            ? fCollectedCharge
            : depositedCharge;

    if (upsetCharge >= fCriticalCharge)
        fUpsetCount += fEventWeight;

    if (fCSV.is_open())
    {
        fCSV
            << event->GetEventID() << ","
            << depositedCharge_fC << ","
            << collectedCharge_fC << ","
            << fTotalEdep / keV << ","
            << fProtonEdep / keV << ","
            << fElectronEdep / keV << ","
            << fRecoilEdep / keV << ","
            << upsetCharge / (1.0e-15 * CLHEP::coulomb) << ","
            << fEventWeight
            << G4endl;
    }

    if (fVerbose)
    {
        G4cout
            << "Event " << event->GetEventID()
            << " | Hits: " << fHits.size()
            << " | Total: " << fTotalEdep / keV << " keV"
            << " | Proton: " << fProtonEdep / keV << " keV"
            << " | Electron: " << fElectronEdep / keV << " keV"
            << " | Recoil: " << fRecoilEdep / keV << " keV"
            << " | Deposited Q: " << depositedCharge_fC << " fC"
            << " | Collected Q: " << collectedCharge_fC << " fC"
            << " | Upset: "
            << (upsetCharge >= fCriticalCharge ? 1 : 0)
            << G4endl;
    }
}

void EventAction::AddHit(const Hit& hit)
{
    fHits.push_back(hit);
}

void EventAction::AddProtonEdep(G4double edep)
{
    fProtonEdep += edep;
    fTotalEdep += edep;
}

void EventAction::AddElectronEdep(G4double edep)
{
    fElectronEdep += edep;
    fTotalEdep += edep;
}

void EventAction::AddRecoilEdep(G4double edep)
{
    fRecoilEdep += edep;
    fTotalEdep += edep;
}

void EventAction::UpdateEventWeight(G4double edep, G4double weight)
{
    if (edep > fMaxSingleEdep)
    {
        fMaxSingleEdep = edep;
        fEventWeight = weight;
    }
}

void EventAction::SetCriticalCharge(G4double qc)
{
    fCriticalCharge = qc * 1.0e-15 * CLHEP::coulomb;
}

void EventAction::SetVerbose(G4bool val)
{
    fVerbose = val;
}

void EventAction::AddCollectedCharge(G4double q)
{
    fCollectedCharge += q;
}

void EventAction::MergeThreadOutputs()
{
    namespace fs = std::filesystem;

    const fs::path dir(kResultsDir);
    const std::string prefix = "events_t";
    const std::string suffix = ".csv";

    std::vector<fs::path> threadFiles;

    for (const auto& entry : fs::directory_iterator(dir))
    {
        const std::string name = entry.path().filename().string();

        if (name.size() > prefix.size() + suffix.size()
            && name.compare(0, prefix.size(), prefix) == 0
            && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
        {
            threadFiles.push_back(entry.path());
        }
    }

    std::ofstream merged(dir / "events.csv");
    bool headerWritten = false;

    for (const auto& file : threadFiles)
    {
        std::ifstream in(file);
        std::string line;
        bool firstLine = true;

        while (std::getline(in, line))
        {
            if (firstLine)
            {
                firstLine = false;

                if (headerWritten)
                    continue; // skip repeated header from later files

                headerWritten = true;
            }

            merged << line << "\n";
        }
    }

    merged.close();

    for (const auto& file : threadFiles)
        fs::remove(file);

    G4cout << "EventAction::MergeThreadOutputs() -- merged "
           << threadFiles.size() << " per-thread file(s) into "
           << (dir / "events.csv").string() << G4endl;
}
