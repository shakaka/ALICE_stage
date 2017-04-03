#ifndef AliAnalysisMuonsCharge_H
#define AliAnalysisMuonsCharge_H
#include "AliAnalysisTaskSE.h"

class TList;
class TObjArray;
class AliMuonEventCuts;
class AliMuonTrackCuts;
class AliCounterCollection;
class AliVParticle;
class AliAODEvent;
class TLorentzVector;

class AliAnalysisMuonsCharge : public AliAnalysisTaskSE {

  public:

    AliAnalysisMuonsCharge();
    AliAnalysisMuonsCharge(const char *name);
    virtual ~AliAnalysisMuonsCharge();

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
      kDimuonPt = 0, ///< dimuon pt
      kMuonPt = 1, ///< single muon pt
      kEventZv = 2,
      kMuonEta = 3,
      kMuonRab = 4,
      kMuonPhi = 5,
      kMuonpDCA = 6,
      kDiMuY = 7,
      kDiMuPhi = 8,
      kDiMuM = 9,
      kDiMuZv = 10,
      kDiMuCh = 11
    };

    AliAODEvent* fAODEvent;       // AOD event

    AliMuonEventCuts *fMuonEventCuts; //< Event cuts
    AliMuonTrackCuts *fMuonTrackCuts; //< Track cuts
    TObjArray *fOutput;               //!< List of histograms for data

    AliCounterCollection *fEventCounters; //!< Event statistics

    ClassDef(AliAnalysisMuonsCharge, 1);

};

#endif
