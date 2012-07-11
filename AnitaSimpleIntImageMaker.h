//////////////////////////////////////////////////////////////////////////////
/////  AnitaSimpleIntImageMaker.h      ANITA Interferrometric Image Maker/////
/////                                                                    /////
/////  Description:                                                      /////
/////     The ANITA class for making interferrometric images             /////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#ifndef ANITASIMPLEINTIMAGEMAKER_H
#define ANITASIMPLEINTIMAGEMAKER_H

//Includes
#include <TObject.h>

#define NUM_PHI_BINS 360
#define NUM_THETA_BINS 180

#define NUM_PAIRS 48

#include "UsefulAnitaEvent.h"

class TH2D;


class AnitaSimpleIntImageMaker : public TObject
{
  AnitaSimpleIntImageMaker(); /// Default constructor
  ~AnitaSimpleIntImageMaker(); // Default destructor


   void fillAntennaPositions();
   void fillAntennaPairs();
   void fillDeltaTArrays();
   Double_t getDeltaTExpected(Int_t ant1, Int_t ant2, Double_t phiWaveRad, Double_t thetaWaveRad);

   TH2D *getInterferometricMap(UsefulAnitaEvent *usefulEvent,AnitaPol::AnitaPol_t pol=AnitaPol::kVertical);

   //Instance generator
   static AnitaSimpleIntImageMaker*  Instance();
   
   static AnitaSimpleIntImageMaker *fgInstance;  
   // protect against multiple instances


   static TGraph* getNormalisedGraph(TGraph *grIn);
   static Double_t fastEvalForEvenSampling(TGraph *grIn, Double_t xvalue);
   static Int_t getPhiSector(Double_t phiWaveDeg);


   Double_t fPhiWaveDeg[NUM_PHI_BINS];
   Double_t fThetaWaveDeg[NUM_PHI_BINS];
   Double_t fPhiWave[NUM_PHI_BINS];
   Double_t fThetaWave[NUM_PHI_BINS];


   Int_t fNumPairs;
   Int_t fFirstAnt[NUM_PAIRS];
   Int_t fSecondAnt[NUM_PAIRS];

   Int_t fPairsPerPhiSector;
   Int_t fPairList[PHI_SECTORS][7];
     

   Double_t fDeltaT[NUM_PAIRS][NUM_PHI_BINS][NUM_THETA_BINS];

   Double_t fZArray[NUM_SEAVEYS];
   Double_t fRArray[NUM_SEAVEYS];
   Double_t fPhiArray[NUM_SEAVEYS];

   ClassDef(AnitaSimpleIntImageMaker,1);
};


#endif //ANITASIMPLEINTIMAGEMAKER_H
