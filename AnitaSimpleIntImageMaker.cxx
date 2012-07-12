#include "AnitaSimpleIntImageMaker.h"
#include "AnitaGeomTool.h"
#include "FFTtools.h"
#include "TH2.h"
#include "TMath.h"
#include <iostream>

ClassImp(AnitaSimpleIntImageMaker);

//Ensure the default pointer is zeroed at start
AnitaSimpleIntImageMaker *AnitaSimpleIntImageMaker::fgInstance=0;


AnitaSimpleIntImageMaker::AnitaSimpleIntImageMaker()
{
  //Default Constructor
  fillAntennaPositions();
  fillAntennaPairs();
  fillDeltaTArrays();
}


AnitaSimpleIntImageMaker::~AnitaSimpleIntImageMaker()
{
  //Default Destructor
}



//______________________________________________________________________________
AnitaSimpleIntImageMaker*  AnitaSimpleIntImageMaker::Instance()
{
  //static function
  if(fgInstance)
    return fgInstance;

  fgInstance = new AnitaSimpleIntImageMaker();
  return fgInstance;
}

void AnitaSimpleIntImageMaker::fillAntennaPositions()
{
  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  for(int ant=0;ant<NUM_SEAVEYS;ant++) {
    fZArray[ant]=fGeomTool->getAntZ(ant);
    fRArray[ant]=fGeomTool->getAntR(ant);
    fPhiArray[ant]=fGeomTool->getAntPhiPositionRelToAftFore(ant);
  }
}


void AnitaSimpleIntImageMaker::fillAntennaPairs()
{
  Int_t pairIndex=0;
  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();

  //Up down ones first
  for(int phi=0;phi<PHI_SECTORS;phi++) {
    int topAnt=fGeomTool->getAntFromPhiRing(phi,AnitaRing::kUpperRing);
    int midAnt=fGeomTool->getAntFromPhiRing(phi,AnitaRing::kLowerRing);
    fFirstAnt[pairIndex]=topAnt;
    fSecondAnt[pairIndex]=midAnt;
    pairIndex++;
  }


  for(int phi=0;phi<PHI_SECTORS;phi++) {
    int rightPhi=(phi+1)%PHI_SECTORS;
    int leftAnt=fGeomTool->getAntFromPhiRing(phi,AnitaRing::kUpperRing);
    int rightAnt=fGeomTool->getAntFromPhiRing(rightPhi,AnitaRing::kUpperRing);
    fFirstAnt[pairIndex]=leftAnt;
    fSecondAnt[pairIndex]=rightAnt;
    pairIndex++;
  }

  for(int phi=0;phi<PHI_SECTORS;phi++) {
    int rightPhi=(phi+1)%PHI_SECTORS;
    int leftAnt=fGeomTool->getAntFromPhiRing(phi,AnitaRing::kLowerRing);
    int rightAnt=fGeomTool->getAntFromPhiRing(rightPhi,AnitaRing::kLowerRing);
    fFirstAnt[pairIndex]=leftAnt;
    fSecondAnt[pairIndex]=rightAnt;
    pairIndex++;
  }

  fNumPairs=pairIndex;
  std::cout << "Setup " << fNumPairs << " antenna pairs out of " << NUM_PAIRS << "\n";

  //Not sure how to make this bit look less crap
  fPairsPerPhiSector=7;
  for(int phi=0;phi<PHI_SECTORS;phi++) {
    int leftPhi=(phi-1)%PHI_SECTORS;
    int rightPhi=(phi+1)%PHI_SECTORS;
    fPairList[phi][0]=phi;
    fPairList[phi][1]=leftPhi;
    fPairList[phi][2]=rightPhi;
    fPairList[phi][3]=leftPhi+16;
    fPairList[phi][4]=phi+16;
    fPairList[phi][5]=leftPhi+32;
    fPairList[phi][6]=phi+32;
  }

}



void AnitaSimpleIntImageMaker::fillDeltaTArrays()
{
  
  Double_t thetaStepDeg=180/NUM_THETA_BINS;
  Double_t phiStepDeg=360/NUM_PHI_BINS;
  
  for(int thetaInd=0;thetaInd<NUM_THETA_BINS;thetaInd++) {
    for(int phiInd=0;phiInd<NUM_PHI_BINS;phiInd++) {
      fPhiWaveDeg[phiInd]=phiInd*phiStepDeg;
      fPhiWave[phiInd]=fPhiWaveDeg[phiInd]*TMath::DegToRad();
      fThetaWaveDeg[thetaInd]=-90+(thetaInd*thetaStepDeg);
      fThetaWave[thetaInd]=fThetaWaveDeg[thetaInd]*TMath::DegToRad();

      for(int pair=0;pair<fNumPairs;pair++) {
	fDeltaT[pair][phiInd][thetaInd]=getDeltaTExpected(fFirstAnt[pair],fSecondAnt[pair],fPhiWave[phiInd],fThetaWave[thetaInd]);
      }
    }
  }
	 
}

Double_t AnitaSimpleIntImageMaker::getDeltaTExpected(Int_t ant1, Int_t ant2, Double_t phiWaveRad, Double_t thetaWaveRad)
{

   Double_t tanThetaW=TMath::Tan(thetaWaveRad);
   
   //   if(tanThetaW>1000) tanThetaW=1000; // Just a little check
   Double_t part1=fZArray[ant1]*tanThetaW - fRArray[ant1] * TMath::Cos(phiWaveRad-fPhiArray[ant1]);
   Double_t part2=fZArray[ant2]*tanThetaW - fRArray[ant2] * TMath::Cos(phiWaveRad-fPhiArray[ant2]);
   
   Double_t tdiff=1e9*((TMath::Cos(thetaWaveRad) * (part1 - part2))/TMath::C());    //returns time in ns
   
   //   std::cerr << tanThetaW << "\t" << tdiff << "\n";
   return tdiff;


}



TH2D *AnitaSimpleIntImageMaker::getInterferometricMap(UsefulAnitaEvent *usefulEvent,AnitaPol::AnitaPol_t pol)
{

  Double_t scale=1./fPairsPerPhiSector;
  TH2D *histMap = new TH2D("histMap","histMap",NUM_PHI_BINS,-180,180,NUM_THETA_BINS,-90,90);

  TGraph *grRaw[NUM_SEAVEYS]={0};
  TGraph *grInt[NUM_SEAVEYS]={0};
  TGraph *grNorm[NUM_SEAVEYS]={0};
  TGraph *grCor[NUM_PAIRS]={0};
  
  for(int ant=0;ant<NUM_SEAVEYS;ant++) {
    //RJN delete this line to include nadirs
    if(ant>=32) continue;
    grRaw[ant]=usefulEvent->getGraph(ant,pol);
    grInt[ant]=FFTtools::getInterpolatedGraph(grRaw[ant],1./2.6);
    grNorm[ant]=getNormalisedGraph(grInt[ant]);
  }


  for(int pair=0;pair<fNumPairs;pair++) {
    int ant1=fFirstAnt[pair];
    int ant2=fSecondAnt[pair];
    grCor[pair]=FFTtools::getCorrelationGraph(grNorm[ant1],grNorm[ant2]);
  }

     
  for(int phiBin=0;phiBin<NUM_PHI_BINS;phiBin++) {
    Int_t whichPhi=getPhiSector(fPhiWaveDeg[phiBin]);
    for(int thetaBin=0;thetaBin<NUM_THETA_BINS;thetaBin++) {
      //I think this is the correct equation to work out the bin number
      //Could just use TH2::GetBin(binx,biny) but the below should be faster
      Int_t globalBin=(phiBin+1)+(thetaBin+1)*(NUM_PHI_BINS+2);
      
      for(int pairInd=0;pairInd<fPairsPerPhiSector;pairInd++) {
	Int_t whichPair=fPairList[whichPhi][pairInd];
	Double_t dt=fDeltaT[whichPair][phiBin][thetaBin];
	Double_t corVal=fastEvalForEvenSampling(grCor[whichPair],dt);
	corVal*=scale;
	Double_t binVal=histMap->GetBinContent(globalBin);
	histMap->SetBinContent(globalBin,binVal+corVal);	
      }
    }
  }
      
  for(int i=0;i<NUM_SEAVEYS;i++) {
    if(grRaw[i]) delete grRaw[i];
    if(grInt[i]) delete grInt[i];
    if(grNorm[i]) delete grNorm[i];
  }
  for(int i=0;i<fNumPairs;i++) {
    if(grCor[i]) delete grCor[i];
  }

  return histMap;
}



TGraph* AnitaSimpleIntImageMaker::getNormalisedGraph(TGraph *grIn)
{
  Double_t rms=grIn->GetRMS(2);
  Double_t mean=grIn->GetMean(2);
  Double_t *xVals = grIn->GetX();
  Double_t *yVals = grIn->GetY();
  Int_t numPoints = grIn->GetN();
  Double_t *newY = new Double_t [numPoints];
  for(int i=0;i<numPoints;i++) {
    newY[i]=(yVals[i]-mean)/rms;
  }
  TGraph *grOut = new TGraph(numPoints,xVals,newY);
  delete [] newY;
  return grOut;
}

Double_t AnitaSimpleIntImageMaker::fastEvalForEvenSampling(TGraph *grIn, Double_t xvalue)
{
  Int_t numPoints=grIn->GetN();
  if(numPoints<2) return 0;
  Double_t *xVals=grIn->GetX();
  Double_t *yVals=grIn->GetY();
  Double_t dx=xVals[1]-xVals[0];
  if(dx<=0) return 0;

  Int_t p0=Int_t((xvalue-xVals[0])/dx);
  if(p0<0) p0=0;
  if(p0>=numPoints) p0=numPoints-2;
  return FFTtools::simpleInterploate(xVals[p0],yVals[p0],xVals[p0+1],yVals[p0+1],xvalue);
}


Int_t AnitaSimpleIntImageMaker::getPhiSector(Double_t phiWaveDeg)
{
  if(phiWaveDeg>(360-11.25) || phiWaveDeg<11.25)
    return 1;
  for(int phi=2;phi<16;phi++) {
    if(phiWaveDeg< 11.25 + (phi-1)*22.5)
      return phi;
  }
  return 0;
}
