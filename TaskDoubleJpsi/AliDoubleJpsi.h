#ifndef AliDoubleJpsi_H
#define AliDoubleJpsi_H
#include "AliAnalysisTaskSE.h"

class TList;
class TObjArray;
class AliMuonEventCuts;
class AliMuonTrackCuts;
class AliCounterCollection;
class AliVParticle;
class AliAODEvent;
class TLorentzVector;

class AliDoubleJpsi : public AliAnalysisTaskSE {

  public:

    AliDoubleJpsi();
    AliDoubleJpsi(const char *name);
    virtual ~AliDoubleJpsi();

    virtual void NotifyRun();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    /// Get track cuts
    AliMuonTrackCuts* GetTrackCuts() { return fMuonTrackCuts; }
    /// Get event cuts
    AliMuonEventCuts* GetEventCuts() { return fMuonEventCuts;}



 private:

    //AliAnalysisTaskSimplePt(const AliAnalysisTaskSimplePt&);
    // AliAnalysisTaskSimplePt& operator=(const AliAnalysisTaskSimplePt&);

    enum eList {
      kMaDbJpsi = 0,
      kPtDbJpsi = 1,
      kYDbJpsi = 2,
      kPhiDbJpsi = 3
    };

    AliAODEvent* fAODEvent;       // AOD event

    AliMuonEventCuts *fMuonEventCuts; //< Event cuts
    AliMuonTrackCuts *fMuonTrackCuts; //< Track cuts
    TObjArray *fOutput;               //!< List of histograms for data

    AliCounterCollection *fEventCounters; //!< Event statistics

    ClassDef(AliDoubleJpsi, 1);

};

#endif
