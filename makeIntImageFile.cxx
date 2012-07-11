#include "AnitaConventions.h"
#include "UsefulAnitaEvent.h"
#include "RawAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "AnitaGeomTool.h"
#include "Adu5Pat.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>

void makeGordonTestFile(int run, int numEnts, char *baseDir, char *outDir=0);

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
  makeGordonTestFile(run,numEnts,"/unix/anita1/flight0809/root/",".");
  stopy.Stop();
  std::cout << "Run " << run << "\n";
  std::cout << "CPU Time: " << stopy.CpuTime() << "\t" << "Real Time: "
	    << stopy.RealTime() << "\n";

}

#define NUM_ANTENNAS 48
#define NUM_SAMPLES 1024
#define NOMINAL_SAMPLING (1./2.6)  //in ns
#define SIGMA_MV 5
  
void makeGordonTestFile(int run, int numEnts, char *baseDir, char *outDir) {
  
  Short_t dataArray[NUM_ANTENNAS][NUM_SAMPLES];

  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  char gpsName[FILENAME_MAX];
  char outName[FILENAME_MAX];
  
  //The locations of the event, header and gps files 
  // The * in the evnt file name is a wildcard for any _X files  
  sprintf(eventName,"%s/run%d/calEventFile%d*.root",baseDir,run,run);
  sprintf(headerName,"%s/run%d/headFile%d.root",baseDir,run,run);
  sprintf(gpsName,"%s/run%d/gpsEvent%d.root",baseDir,run,run);

  Int_t useCalibratedFiles=0;

  //Define and zero the class pointers
  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
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

  TFile *fpOut = new TFile("temp.root","RECREATE");


  Long64_t maxEntry=headTree->GetEntries(); 
  if(numEnts && maxEntry>numEnts) maxEntry=numEnts;

  Int_t starEvery=maxEntry/100;
  if(starEvery==0) starEvery=1;
  
  std::cout <<  "There are " << maxEntry << " events to proces\n";
  Long64_t countEvents=0;
  for(Long64_t entry=0;entry<maxEntry;entry++) {
     if(entry%starEvery==0) std::cerr << "*";

     //Get header
     headTree->GetEntry(entry);

     //    if( (header->triggerTimeNs>0.4e6) || (header->triggerTimeNs<0.25e6) )  
     //Now cut to only process the Taylor Dome pulses
     if( (header->triggerTimeNs>3e6 || header->triggerTimeNs<1e5) )
       continue; 

     
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
     memset(dataArray,0,sizeof(Short_t)*NUM_SAMPLES*NUM_ANTENNAS);

     for(int ring=AnitaRing::kUpperRing;ring<AnitaRing::kNotARing;ring++) {
       for(int phi=0;phi<16;phi++) {
	  Int_t chanIndex=AnitaGeomTool::getChanIndexFromRingPhiPol(AnitaRing::AnitaRing_t(ring),phi,AnitaPol::kVertical);
	 TGraph *gr = realEvent->getGraph((AnitaRing::AnitaRing_t)ring,phi,AnitaPol::kVertical);
	 if(gr) {
	    if(ring==0 && phi==0)
	    Int_t numPoints=gr->GetN();
	    if(countEvents==0) {
	       Double_t maxVal=0;
	       Double_t maxTime=0;
	       for(int samp=0;samp<numPoints;samp++) {
		  if(gr->GetY()[samp]>maxVal) {
		     maxVal=gr->GetY()[samp];
		     maxTime=gr->GetX()[samp];
		  }
	       }
	       std::cout << ring << "\t" << phi << "\t" << chanIndex << "\t" << maxVal << "\t" << maxTime << "\n";
	    }
	   //	   std::cout << ring << "\t" << phi << "\t" << chanIndex <<  "\n";
	   Double_t startTime= gr->GetX()[0];
	   Int_t startInd=TMath::Nint(startTime/NOMINAL_SAMPLING)+512;
	   //       std::cout << startTime << "\t" << startInd << "\n";
	   for(int samp=0;samp<NUM_SAMPLES;samp++) {
	     if(samp<startInd || samp>=startInd+numPoints)
	       dataArray[16*ring+phi][samp]=TMath::Nint(gRandom->Gaus(0,SIGMA_MV));
	     else 
	       dataArray[16*ring+phi][samp]=TMath::Nint(gr->GetY()[samp-startInd]);
	   }
	   delete gr;
	 }	 
	 else {
	   for(int samp=0;samp<NUM_SAMPLES;samp++) {
	     dataArray[16*ring+phi][samp]=TMath::Nint(gRandom->Gaus(0,SIGMA_MV));
	   }
	   
	 }
       } 
	 
     }
     
     //     std::cout << run << "\t" << countEvents << "\t" << header->eventNumber << "\t" << realEvent->eventNumber << "\n";

     sprintf(outName,"/unix/anita2/triggerTest/dataFile_%d_%d.txt",run,countEvents);
     std::ofstream Output(outName);
     for(int ant=0;ant<NUM_ANTENNAS;ant++) {
       for(int samp=0;samp<NUM_SAMPLES;samp++) {
	 Output << dataArray[ant][samp] << "\t";
       }
       Output << "\n";
     }
     Output.close();
     //     return ;
  

     if(realEvent) delete realEvent;

     countEvents++;
  }
  std::cerr << "\n";
  std::cout << countEvents << " events processed\n";
}
     
