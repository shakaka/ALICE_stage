// For more information on the TSelector framework see 
// $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The file for this selector can be found at
// http://lcg-heppkg.web.cern.ch/lcg-heppkg/ROOT/eventdata.root
// i.e run
//   root [0] f = TFile::Open("http://lcg-heppkg.web.cern.ch/lcg-heppkg/ROOT/eventdata.root");
//   root [1] EventTree->Process("EventSelector.C+")

// The following methods are defined in this file:
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers, a convenient place to create your histograms.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("EventSelector.C")
// Root > T->Process("EventSelector.C","some options")
// Root > T->Process("EventSelector.C+")
//

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"

const Int_t kMaxfParticles = 1293;

class EventSelector : public TSelector {
public :
   // Variables used to store the data
   Int_t       fNumberOfEvents; // Total number of events
  // Variables used to access and store the data
   TTreeReader fReader;                       // The tree reader 
   TTreeReaderValue<Int_t> fCurrentEventSize; // Size of the current event

   EventSelector(TTree * = 0):
      fTotalDataSize(0),
      fCurrentEventSize(fReader, "fEventSize") { }
   virtual ~EventSelector() { }

   virtual void    Init(TTree *tree);


   // EventSelector(TTree * = 0) { }
   // virtual ~EventSelector() { }

   // virtual void    Init(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual Bool_t  Process(Long64_t entry);
   virtual void    Terminate();
   virtual Int_t   Version() const { return 2; }

   Int_t fTotalDataSize;   // Sum of data size (in bytes) of all events
   


   hPosX = new TH1F("hPosX", "Position in X", 20, -5, 5);


   ClassDef(EventSelector,0);
};

void EventSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
  fReader.SetTree(tree);
}

void EventSelector::SlaveBegin(TTree *tree)
{
   // SlaveBegin() is a good place to create histograms. 
   // For PROOF, this is called for each worker.
   // The TTree* is there for backward compatibility; e.g. PROOF passes 0.

}

Bool_t EventSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree to be 
   // processed. The entry argument specifies which entry in the currently
   // loaded tree is to be processed.
   // It can be passed to either EventSelector::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the TTree.
   //
   // This function should contain the "body" of the analysis: select relevant
   // tree entries, run algorithms on the tree entry and typically fill histograms.
  

   // print some information about the current entry
   printf("Processing Entry number %lld\n", entry);
   // increase the total number of entries
   ++fNumberOfEvents;


   // Tell the TTree reader to get the data for
   // the entry number "entry" in the current tree:
   fReader.SetLocalEntry(entry);

   // We can still print some informations about the current event
   //printf("Size of Event %ld = %d Bytes\n", entry, *fCurrentEventSize);

   // compute the total size of all events; dereference the TTreeReaderValue
   // using '*' to get the value it refers to, just like an iterator.
   fTotalDataSize += *fCurrentEventSize;

   hPosX->Sumw2();


   while (fReader.Next()){
     for(int i=0; i<particlesMomentum.GetSize(); i++){
        if (particlesMomentum[i] > 40.0)
	  hPosX->Fill(particlesPosX[i]);
     }
   }
   hPosX->Fit("pol2");
   hPosX->Draw();






   return kTRUE;
}

void EventSelector::Terminate()
{
   // The Terminate() function is the last function to be called during the
   // analysis of a tree with a selector. It always runs on the client, it can
   // be used to present the results graphically or save the results to file.


   // print the result
   printf("\nTotal Number of Events: %d\n", fNumberOfEvents);

   int sizeInMB = fTotalDataSize/1024/1024;
   printf("Total size of all events: %d MB\n", sizeInMB);
}
