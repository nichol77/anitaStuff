#include "AnitaConventions.h"
#include "UsefulAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "AnitaGeomTool.h"
#include "AnitaSimpleIntImageMaker.h"
#include "Adu5Pat.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>

void makeIntImageFile(int run, int numEnts, const char *baseDir, const char *outDir=0);

int main(int argc, char **argv) {
  int run=165;
  if(argc>1) {
    run=atoi(argv[1]);
  }
  int numEnts=0;
  if(argc>2) {
    numEnts=atoi(argv[2]);
  }
  std::cout << "Making correlation summary tree for run: " << run << "\n";

  TStopwatch stopy;
  stopy.Start();
  makeIntImageFile(run,numEnts,"/unix/anita1/flight0809/root/",".");
  stopy.Stop();
  std::cout << "Run " << run << "\n";
  std::cout << "CPU Time: " << stopy.CpuTime() << "\t" << "Real Time: "
	    << stopy.RealTime() << "\n";

}

void makeIntImageFile(int run, int numEnts, const char *baseDir, const char *outDir) {

  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  char histName[FILENAME_MAX];
  
  //The locations of the event, header and gps files 
  // The * in the evnt file name is a wildcard for any _X files  
  sprintf(eventName,"%s/run%d/calEventFile%d*.root",baseDir,run,run);
  sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);
  sprintf(gpsName,"%s/run%d/gpsEvent%d.root",baseDir,run,run);
  sprintf(outName,"%s/mapFile%d.root",outDir,run);
  Int_t useCalibratedFiles=0;

  //Define and zero the class pointers
  //  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  RawAnitaEvent *event = 0;
  CalibratedAnitaEvent *calEvent = 0;
  RawAnitaHeader *header =0;
  
  //Need a TChain for the event files as they are split in to sub files to ensure no single file is over 2GB
  TChain *eventChain = new TChain("eventTree");
  eventChain->Add(eventName);

  //Here we check if there are any entries in the tree of CalibratedAnitaEvent objects
  if(eventChain->GetEntries()>0) {
    //If there are entries we can set the branc address
    eventChain->SetBranchAddress("event",&calEvent);
    useCalibratedFiles=1;
    std::cout << "Using calibrated event files\n";
  }
  else {
    //If there aren't any entries we can try raw event files instead
    sprintf(eventName,"%s/run%d/eventFile%d*.root",baseDir,run,run);
    eventChain->Add(eventName);
    eventChain->SetBranchAddress("event",&event);
    std::cout << "Using raw event files\n";
  }
  
  //Now open the header and GPS files
  TFile *fpHead = TFile::Open(headerName);
  TTree *headTree = (TTree*) fpHead->Get("headTree");
  headTree->SetBranchAddress("header",&header);

  TFile *fpOut = new TFile(outName,"RECREATE");
  AnitaSimpleIntImageMaker *imageMaker = AnitaSimpleIntImageMaker::Instance();

  Long64_t maxEntry=headTree->GetEntries(); 
  if(numEnts && maxEntry>numEnts) maxEntry=numEnts;

  Int_t starEvery=maxEntry/100;
  if(starEvery==0) starEvery=1;
  
  std::cout <<  "There are " << maxEntry << " events to proces\n";
  for(Long64_t entry=0;entry<maxEntry;entry++) {
     if(entry%starEvery==0) std::cerr << "*";

     //Get header
     headTree->GetEntry(entry);

     //Now cut to only process the Taylor Dome pulses
     //     if( (header->triggerTimeNs>3e6 || header->triggerTimeNs<1e5) )
     //       continue; 

     
     //Get event
     eventChain->GetEntry(entry);

     //Now we can make a UsefulAnitaEvent (or a UsefulAnitaEvent)
     UsefulAnitaEvent *realEvent=0;
     if(useCalibratedFiles) {
       //If we have CalibratedAnitaEvent then the constructor is just
       realEvent = new UsefulAnitaEvent(calEvent);
     }
     else {
       //If we have RawAnitaEvent then we have to specify the calibration option
       realEvent = new UsefulAnitaEvent(event,WaveCalType::kDefault,header);
     }
     TH2D *map = imageMaker->getInterferometricMap(realEvent,AnitaPol::kVertical);
     sprintf(histName,"histMapV_%d",realEvent->eventNumber);
     map->SetNameTitle(histName,histName);
     //     map->Write();
     delete map;
  }
  fpOut->Close();
}
     
