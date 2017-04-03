#ifndef AliAnalysisJpsi_H
#define AliAnalysisJpsi_H
#include "AliAnalysisTaskSE.h"

class TList;
class TObjArray;
class AliMuonEventCuts;
class AliMuonTrackCuts;
class AliCounterCollection;
class AliVParticle;
class AliAODEvent;
class TLorentzVector;

class AliAnalysisJpsi : public AliAnalysisTaskSE {

  public:

    AliAnalysisJpsi();
    AliAnalysisJpsi(const char *name);
    virtual ~AliAnalysisJpsi();

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
      kDiMuPt = 0, ///< dimuon pt
      kDiMuY = 1,
      kDiMuPhi = 2,
      kDiMuM = 3,
      kDiMuZv = 4,
      kDiMuCh = 5
    };

    AliAODEvent* fAODEvent;       // AOD event

    AliMuonEventCuts *fMuonEventCuts; //< Event cuts
    AliMuonTrackCuts *fMuonTrackCuts; //< Track cuts
    TObjArray *fOutput;               //!< List of histograms for data

    AliCounterCollection *fEventCounters; //!< Event statistics

    ClassDef(AliAnalysisJpsi, 1);

};

#endif
