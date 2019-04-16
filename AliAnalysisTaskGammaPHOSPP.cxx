#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TList.h>
#include <THashList.h>
#include <TClonesArray.h>  
#include <TGeoGlobalMagField.h>
#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TClonesArray.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskGammaPHOSPP.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODTrack.h"
#include "AliPHOSAodCluster.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliTriggerAnalysis.h"
#include "AliVVZERO.h"
#include "AliPHOSTriggerUtils.h"
#include "AliAODMCParticle.h"
#include "AliPHOSTenderSupply.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliInputEventHandler.h"
#include "AliTriggerAnalysis.h"
#include "AliOADBContainer.h"

#include "AliAnalysisManager.h"

// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Yuri Kharlov
// Date   : 28.05.2009

ClassImp(AliAnalysisTaskGammaPHOSPP)

//________________________________________________________________________
AliAnalysisTaskGammaPHOSPP::AliAnalysisTaskGammaPHOSPP(const char *name) 
: AliAnalysisTaskSE(name),
  fESDtrackCuts(0),
  fOutputContainer(0),
  fOutputContainer2(0),
  fPHOSEvent(0),
  fnCINT1B(0),
  fnCINT1A(0),
  fnCINT1C(0),
  fnCINT1E(0),
  fPHOSGeo(0),
  fBCgap(525e-09),
  fEventCounter(0),
  fAllEventCounter(0),
  fEventCentrality(0),
  fInPHOS(0),
  fInPHOS2(0),
  fEvent(0),
  fMCArray(0),
  fPIDResponse(0x0), //!
  fWeightFunction(0),
  fWeightFunction2(0),
  fCurrFileName(0), 
  fCheckMCCrossSection(kFALSE),
  fh1Xsec(0),      
  fh1Trials(0),
  fAvgTrials(-1),
  fTriggerAnalysis(new AliTriggerAnalysis)
{
  
  // Output slots #0 write into a TH1 container
  DefineOutput(1,THashList::Class());
  DefineOutput(2,THashList::Class());

// Constructor
  Int_t nBin=10 ;
  for(Int_t i=0; i<nBin; i++)
  {
    for(Int_t j=0; j<2; j++)
         fPHOSEvents[i][j] = 0 ;
   }
  
  Double_t fVtx0[3], fVtx5[3];   
  for(Int_t i = 0; i < 3; i++)
  {
    fVtx0[i] = 0;
    fVtx5[i] = 0;
  } 
    
   fWeightFunction2 = new TF1("fWeightFunction2", "([0] + [1]*x + [2]*x*x)/(1. + [3]*x + [4]*x*x) + [5]*x", 0.1, 400); 
}

//________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::UserCreateOutputObjects()
{
  // Create histograms, called once

  // AOD histograms
  if(fOutputContainer != NULL)
  {
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  if(fOutputContainer2 != NULL)
  {
    delete fOutputContainer2;
  }
  fOutputContainer2 = new THashList();
  fOutputContainer2->SetOwner(kTRUE);
  

  const Int_t Ncut = 4;
  TString cut[Ncut] = {"all", "cpv", "disp", "both"};

 //PID Response 
     AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
     if(!man)
         AliFatal("Could not find manager");
     AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> (man->GetInputEventHandler());
     if(!inputHandler)
         AliFatal("No input event handler");
     fPIDResponse = dynamic_cast<AliPIDResponse *>(inputHandler->GetPIDResponse());
     if (!fPIDResponse)
         AliFatal("PIDResponse object was not created"); // Escalated to fatal. This task is unusable without PID response.
    fPIDResponse -> SetUseTPCMultiplicityCorrection(kFALSE);
    fPIDResponse->SetCurrentMCEvent(MCEvent()); //!!

    fh1Xsec = new TH1F("hXsec","xsec from pyxsec.root",1,0,1);
    fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
    fOutputContainer->Add(fh1Xsec);
    
    fh1Trials = new TH1F("hTrials","trials root file",1,0,1);
    fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
    fOutputContainer->Add(fh1Trials);

  fOutputContainer->Add(new TH1I("hCellMultEvent"  ,"PHOS cell multiplicity per event"    ,2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM1","PHOS cell multiplicity per event, M1",2000,0,2000));
  for(Int_t imod = 1; imod < 5; imod ++)
  {
   fOutputContainer->Add(new TH1I(Form("hCellMultEventM%d", imod),Form("PHOS cell multiplicity per event, M%id", imod), 2000, 0, 2000));
  }

  fOutputContainer->Add(new TH1I("hClusterMult"      ,"CaloCluster multiplicity"     ,100,0,100));
  fOutputContainer->Add(new TH1I("hPHOSClusterMult"  ,"PHOS cluster multiplicity"    ,100,0,100));

  for(Int_t imod = 1; imod < 5; imod ++)
  {
  fOutputContainer->Add(new TH1I(Form("hPHOSClusterMultM%d",imod),Form("PHOS cluster multiplicity, M%d",imod),100,0,100));
  }
  fOutputContainer->Add(new TH1F("hCellEnergy"  ,"Cell energy"            ,5000,0.,50.));
  for(Int_t imod = 1; imod < 5; imod ++)
  {
     fOutputContainer->Add(new TH1F(Form("hCellEnergyM%d",imod),Form("Cell energy in module %d", imod),5000,0.,50.));
  }
 
  fOutputContainer->Add(new TH1F("hClusterEnergy"  ,"Cluster energy"      ,5000,0.,50.));
  for(Int_t imod = 1; imod < 5; imod ++)
  {  
  fOutputContainer->Add(new TH1F(Form("hClusterEnergyM%d", imod),Form("Cluster energy, M%d",imod)  ,5000,0.,50.));
  }


 for(Int_t iCut = 0; iCut < Ncut; iCut++) 
 {
  fOutputContainer->Add(new TH2F(Form("hClusterEvsN_%s", cut[iCut].Data())  ,"Cluster energy vs digit multiplicity"    ,5000,0.,50.,40,0.,40.));
  for(Int_t imod = 1; imod < 5; imod ++)
  {
   fOutputContainer->Add(new TH2F(Form("hClusterEvsN_%s_M%d", cut[iCut].Data(), imod),Form("Cluster energy vs digit multiplicity, M%d",imod),5000,0.,50.,40,0.,40.));
  } 


  fOutputContainer->Add(new TH1I(Form("hCellMultClu_%s", cut[iCut].Data())  ,"Cell multiplicity per cluster"    ,200,0,200));
  for(Int_t imod = 1; imod < 5; imod ++)
  {
  fOutputContainer->Add(new TH1I(Form("hCellMultClu_%s_M%d", cut[iCut].Data(), imod),Form("Cell multiplicity per cluster, M%d", imod),200,0,200));
  }

 }

  fOutputContainer->Add(new TH1I("hModule","Module events", 5, 0., 5.));
  fOutputContainer->Add(new TH1F("hSelEvents","Selected events",9, -0.5, 8.5));
  fOutputContainer->Add(new TH1F("hSelEventspi0","Selected events for pi0", 9, -0.5, 8.5)); 

  for(Int_t imod = 1; imod < 5; imod ++)
  {
  fOutputContainer->Add(new TH2F(Form("hCellNXZM%d", imod), Form("Cell (X,Z), M%d" , imod), 64, 0.5, 64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F(Form("hCellEXZM%d", imod), Form("Cell E(X,Z), M%d", imod), 64, 0.5, 64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F(Form("hCluNXZM%d",  imod), Form("Clu (X,Z), M%d",   imod), 64, 0.5, 64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F(Form("hCluEXZM%d",  imod), Form("Clu E(X,Z), M%d",  imod), 64, 0.5, 64.5, 56,0.5,56.5));
  }

  Int_t nM       = 750;
  Double_t mMin  = 0.0;
  Double_t mMax  = 1.5;
  Int_t nPt      = 400;
  Double_t ptMin = 0;
  Double_t ptMax = 40;
  fOutputContainer->Add(new TH2F("hAsymPtPi0","(A,p_{T})_{#gamma#gamma} #pi^{0}", 20,0., 1.,  40,0., 20.));
  fOutputContainer->Add(new TH2F("hAsymPtEta","(A,p_{T})_{#gamma#gamma} #eta",    20,0., 1.,  40,0., 20.));

  for(Int_t imod = 1; imod < 5; imod ++)
  {
  fOutputContainer->Add(new TH2F(Form("hAsymPtPi0M%d",imod) ,Form("(A,p_{T})_{#gamma#gamma} #pi^{0}. M%d",imod),  20, 0., 1., 40, 0., 20.));
   for(Int_t imodd = imod+1; imodd <5; imodd++)
   {
      fOutputContainer->Add(new TH2F(Form("hAsymPtPi0M%d%d", imod, imodd),Form("(A,p_{T})_{#gamma#gamma} #pi^{0}. M%d%d",imod, imodd) , 20, 0., 1., 40, 0., 20.));
   }
  }


  fOutputContainer->Add(new TH2F("hMassPtA10" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0"   ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA08" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.8"   ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7"   ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA01" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.1"   ,nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, no vtx" ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, no vtx" ,nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, vtx"     ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, vtx"     ,nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, V0AND",nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, V0AND",nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, pileup",nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, pileup",nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMassPtvA10","(M,p_{T})_{#gamma#gamma}, 0<A<1.0, ESD vertex",nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPtvA07","(M,p_{T})_{#gamma#gamma}, 0<A<0.7, ESD vertex",nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH3F("hMassPtCA10","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA10_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA10_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA10_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07","(M,p_{T},C)_{#gamma#gamma}, 0<A<0.7" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));

  fOutputContainer->Add(new TH2F("hMassSingle_all","(M,p_{T})_{#gamma#gamma}, no PID" ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassSingle_cpv","(M,p_{T})_{#gamma#gamma}, no PID" ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassSingle_disp","(M,p_{T})_{#gamma#gamma}, no PID" ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassSingle_both","(M,p_{T})_{#gamma#gamma}, no PID" ,nM, mMin, mMax, nPt, ptMin, ptMax));

  for(Int_t imod = 1; imod < 5; imod ++)
  {  
  fOutputContainer->Add(new TH2F(Form("hMassPtM%d", imod), Form("(M,p_{T})_{#gamma#gamma}, module %d",  imod), nM, mMin, mMax, nPt, ptMin, ptMax));
   for(Int_t imodd = imod+1; imodd < 5; imodd++)
   {
    fOutputContainer->Add(new TH2F(Form("hMassPtM%d%d",imod, imodd),Form("(M,p_{T})_{#gamma#gamma}, modules %d, %d", imod, imodd), nM, mMin, mMax, nPt, ptMin, ptMax));
   }
  }


  fOutputContainer->Add(new TH2F("hMassPt20cm","(M,p_{T})_{#gamma#gamma}, |z|<20 cm"   ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPt40cm","(M,p_{T})_{#gamma#gamma}, 20<|z|<40 cm",nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPt60cm","(M,p_{T})_{#gamma#gamma}, |z|>40 cm"   ,nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMassPtN3","(M,p_{T})_{#gamma#gamma}, N_{cell}>2"  ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPtN4","(M,p_{T})_{#gamma#gamma}, N_{cell}>3"  ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPtN5","(M,p_{T})_{#gamma#gamma}, N_{cell}>4"  ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMassPtN6","(M,p_{T})_{#gamma#gamma}, N_{cell}>5"  ,nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMiAsymPt","(A,p_{T})_{#gamma#gamma}"                ,50,0.,1.,    nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0"   ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA08" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.8"   ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7"   ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA01" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.1"   ,nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, no vtx" ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, no vtx" ,nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, vtx"     ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, vtx"     ,nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, V0AND",nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, V0AND",nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, pileup",nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, pileup",nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtvA10","(M,p_{T})_{#gamma#gamma}, 0<A<1.0, ESD vertex",nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtvA07","(M,p_{T})_{#gamma#gamma}, 0<A<0.7, ESD vertex",nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH3F("hMiMassPtCA10","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA10_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA10_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA10_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07","(M,p_{T},C)_{#gamma#gamma}, 0<A<0.7" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM, mMin, mMax, nPt, ptMin, ptMax,8,0.,8.));

  fOutputContainer->Add(new TH2F("hMiMassSingle_all","(M,p_{T})_{#gamma#gamma}, no PID" ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassSingle_cpv","(M,p_{T})_{#gamma#gamma}, no PID" ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassSingle_disp","(M,p_{T})_{#gamma#gamma}, no PID" ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassSingle_both","(M,p_{T})_{#gamma#gamma}, no PID" ,nM, mMin, mMax, nPt, ptMin, ptMax));

  for(Int_t imod = 1; imod < 5; imod ++)
  {
  fOutputContainer->Add(new TH2F(Form("hMiMassPtM%d", imod),Form("(M,p_{T})_{#gamma#gamma}, module %d",imod) ,nM, mMin, mMax, nPt, ptMin, ptMax));
   for(Int_t imodd = imod+1; imodd <5; imodd ++)
   {
    fOutputContainer->Add(new TH2F(Form("hMiMassPtM%d%d", imod, imodd), Form("(M,p_{T})_{#gamma#gamma}, modules %d, %d", imod, imodd), nM, mMin, mMax, nPt, ptMin, ptMax));
   }
  }


  fOutputContainer->Add(new TH2F("hMiMassPt20cm","(M,p_{T})_{#gamma#gamma}, |z|<20 cm"   ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPt40cm","(M,p_{T})_{#gamma#gamma}, 20<|z|<40 cm",nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPt60cm","(M,p_{T})_{#gamma#gamma}, |z|>40 cm"   ,nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtN3","(M,p_{T})_{#gamma#gamma}, N_{cell}>2"  ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtN4","(M,p_{T})_{#gamma#gamma}, N_{cell}>3"  ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtN5","(M,p_{T})_{#gamma#gamma}, N_{cell}>4"  ,nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtN6","(M,p_{T})_{#gamma#gamma}, N_{cell}>5"  ,nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer->Add(new TH1F("hPhotonKappa","#kappa(#gamma)",400,0.,40.));
  fOutputContainer->Add(new TH1F("hPhotonPt","p_{T}(#gamma)",400,0.,40.));
  fOutputContainer->Add(new TH1F("hPhotonPx","p_{x}(#gamma)",400,0.,40.));
  fOutputContainer->Add(new TH1F("hPhotonPy","p_{y}(#gamma)",400,0.,40.));

  Sumw2Histogram("hPhotonPt");

  fOutputContainer->Add(new TH1F("hTrigClass","Trigger class",5,0.5,5.5));
  fOutputContainer->Add(new TH1F("hNPileupVtx","Number of SPD pileup vertices",10,0.,10.));
  fOutputContainer->Add(new TH1F("hZPileupVtx", "Z pileup", 200, -50., 50.));

  fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200,-50.,+50.));
  fOutputContainer->Add(new TH1F("hNvertexTracks","N of primary tracks from the primary vertex",150,0.,150.));
  fOutputContainer->Add(new TH1F("hTrackMult","Charged track multiplicity",150,0.,150.));

  fOutputContainer->Add(new TH1F("hV0Atime","V0A time",1200,-6.e-6,+6.e-6));
  fOutputContainer->Add(new TH1F("hV0Ctime","V0C time",1200,-6.e-6,+6.e-6));
  fOutputContainer->Add(new TH2F("hV0AV0Ctime","V0A time vs V0C time",120,-6.e-6,+6.e-6 ,120,-6.e-6,+6.e-6));


 for(Int_t iCut = 0; iCut < Ncut; iCut++) 
 {
  fOutputContainer->Add(new TH1F(Form("hClustPt_%s", cut[iCut].Data()),"Cluster P_{t} spectrum",nPt,ptMin,ptMax));
  Sumw2Histogram(Form("hClustPt_%s", cut[iCut].Data()));

  for(Int_t imod = 1; imod <5; imod++)
  {
  fOutputContainer->Add(new TH1F(Form("hClustPt_%s_mod%d", cut[iCut].Data(), imod), Form("Cluster P_{t} spectrum, M%d", imod), nPt, ptMin, ptMax));
  Sumw2Histogram(Form("hClustPt_%s_mod%d", cut[iCut].Data(), imod));  
  }
  
  fOutputContainer->Add(new TH2F(Form("hClustPtvsNcl_%s", cut[iCut].Data()),"Minumum number of cells in a cluster vs cluster P_{t}",nPt,ptMin,ptMax, 5, 3, 8));
  Sumw2Histogram(Form("hClustPtvsNcl_%s", cut[iCut].Data()));


  fOutputContainer->Add(new TH2F(Form("hCentralityvsClustPt_%s", cut[iCut].Data()),"Centrality vs ClustPt", nPt,ptMin,ptMax, 8,0., 8.));
    Sumw2Histogram(Form("hCentralityvsClustPt_%s", cut[iCut].Data()));

 }

  fOutputContainer->Add(new TH1F("hTOF","Time of Flight", 1000, 0, 1e-7));
  fOutputContainer->Add(new TH1I("hBC", "BC", 1000, 0, 1000));

  fOutputContainer->Add(new TH1F("hvt0vsvt5", "Momenta from vertices", 200, -1., 1.));
 
  fOutputContainer2->Add(new TH1F("hWeights","Particle weights", 2000, 0., 2.));

  fOutputContainer->Add(new TH2F("hClustM02","Cluster M02 vs p_{T}",nPt,ptMin,ptMax, 100,0,10));
  fOutputContainer->Add(new TH2F("hClustM20","Cluster M20 vs p_{T}",nPt,ptMin,ptMax, 100,0,10));


  fOutputContainer->Add(new TH2F("hMinvClustPt_all",   "Invarinant mass vs cluster p_{T}",nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMinvClustPt_cpv",   "Invarinant mass vs cluster p_{T}, CPV cut",nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMinvClustPt_disp",  "Invarinant mass vs cluster p_{T}",nM, mMin, mMax, nPt, ptMin, ptMax));
  fOutputContainer->Add(new TH2F("hMinvClustPt_both",  "Invarinant mass vs cluster p_{T}",nM, mMin, mMax, nPt, ptMin, ptMax));

//--------------------Invariant masses, dispersion cut=======

  
   fOutputContainer->Add(new TH2F("hMinv_all_all","Invariant mass vs pt all-all",nM, mMin, mMax, nPt, ptMin, ptMax));
   fOutputContainer->Add(new TH2F("hMinv_all_all_mix","Invariant mass vs pt all-all",nM, mMin, mMax, nPt, ptMin, ptMax));

   fOutputContainer->Add(new TH2F("hMinv_antidisp_all","Invariant mass vs pt antiDisp-all",nM, mMin, mMax, nPt, ptMin, ptMax));
   fOutputContainer->Add(new TH2F("hMinv_antidisp_all_mix","Invariant mass vs pt antiDisp-all",nM, mMin, mMax, nPt, ptMin, ptMax));

   fOutputContainer->Add(new TH2F("hMinv_disp_all","Invariant mass vs pt Disp-all",nM, mMin, mMax, nPt, ptMin, ptMax));
   fOutputContainer->Add(new TH2F("hMinv_disp_all_mix","Invariant mass vs pt Disp-all",nM, mMin, mMax, nPt, ptMin, ptMax));


   fOutputContainer->Add(new TH2F("hMinv_anticpv_all","Invariant mass vs pt antiCPV-all",nM, mMin, mMax, nPt, ptMin, ptMax));
   fOutputContainer->Add(new TH2F("hMinv_anticpv_all_mix","Invariant mass vs pt antiCPV-all",nM, mMin, mMax, nPt, ptMin, ptMax));

   fOutputContainer->Add(new TH2F("hMinv_cpv_all","Invariant mass vs pt CPV-all",nM, mMin, mMax, nPt, ptMin, ptMax));
   fOutputContainer->Add(new TH2F("hMinv_cpv_all_mix","Invariant mass vs pt CPV-all",nM, mMin, mMax, nPt, ptMin, ptMax));

  fOutputContainer ->Add(new TH1F("hTrackCharge","Charge of track", 9, -4.5, 4.5));
  fOutputContainer2->Add(new TH1F("hMCTrackCharge","Charge of MC track", 9, -4.5, 4.5));
  fOutputContainer->Add(new TH2F("hTracksPt2D","Tracks ID vs p_{T}",nPt,ptMin,ptMax, 9,0.,9.));
  fOutputContainer->Add(new TH1F("hTracksPt_beta","Electron tracks p_{T}",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_mu","#mu tracks p_{T}",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_pi","#pi^{#pm} tracks p_{T}",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_ka","K^{#pm} tracks p_{T}",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_p","p, #bar{p} tracks vs p_{T}",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_de","^{2}_{1}H tracks vs p_{T}",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_tri","^{3}_{1}H tracks vs p_{T}",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_he3","^{3}_{2}He tracks vs p_{T}",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hTracksPt_alpha","^{4}_{2} He tracks vs p_{T}",nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hTracksofClusts","Tracks vc p_{T}.",nPt,ptMin,ptMax, 9,0.,9.));
  fOutputContainer->Add(new TH2F("hTracksofClusts_cpv","Tracks vc p_{T}, CPV cut.",nPt,ptMin,ptMax, 9,0.,9.));
  fOutputContainer->Add(new TH2F("hTracksofClusts_disp","Tracks vc p_{T}, shape cut.", nPt,ptMin,ptMax, 9,0.,9.));
  fOutputContainer->Add(new TH2F("hTracksofClusts_both","Tracks vc p_{T}, CPV + shape cuts." ,nPt,ptMin,ptMax, 9,0.,9.));

  
  fOutputContainer->Add(new TH2F("hTracksOfPi_OneSigma", "#pi^{+} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfPi_OneSigma_disp", "#pi^{+} clusters, Disp cut", nPt,ptMin,ptMax,1000.,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPi_OneSigma", "#pi^{-} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPi_OneSigma_disp", "#pi^{-} clusters, Disp cut", nPt,ptMin,ptMax,1000.,0.,100.));

  fOutputContainer->Add(new TH2F("hTracksOfPi_ThreeSigma", "#pi^{+} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfPi_ThreeSigma_disp", "#pi^{+} clusters, Disp cut", nPt,ptMin,ptMax,1000.,0.,100.));  
  fOutputContainer->Add(new TH2F("hTracksOfAntiPi_ThreeSigma", "#pi^{-} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPi_ThreeSigma_disp", "#pi^{-} clusters, Disp cut", nPt,ptMin,ptMax,1000.,0.,100.));

  fOutputContainer->Add(new TH2F("hTracksOfPr_OneSigma", "p clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfPr_OneSigma_disp", "p clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPr_OneSigma", "#bar{p} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPr_OneSigma_disp", "#bar{p} clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));

  fOutputContainer->Add(new TH2F("hTracksOfPr_ThreeSigma", "p clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfPr_ThreeSigma_disp", "p clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPr_ThreeSigma", "#bar{p} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiPr_ThreeSigma_disp", "#bar{p} clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));
 
  fOutputContainer->Add(new TH2F("hTracksOfKa_OneSigma", "K^{+} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfKa_OneSigma_disp", "K^{+} clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiKa_OneSigma", "K^{-} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiKa_OneSigma_disp", "K^{-} clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));

  fOutputContainer->Add(new TH2F("hTracksOfKa_ThreeSigma", "K^{+} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfKa_ThreeSigma_disp", "K^{+} clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiKa_ThreeSigma", "K^{-} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiKa_ThreeSigma_disp", "K^{-} clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));

  fOutputContainer->Add(new TH2F("hTracksOfBeta_OneSigma", "#beta^{-} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfBeta_OneSigma_disp", "#beta^{-} clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiBeta_OneSigma", "#beta^{+} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiBeta_OneSigma_disp", "#beta^{+} clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));

  fOutputContainer->Add(new TH2F("hTracksOfBeta_ThreeSigma", "#beta^{-} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfBeta_ThreeSigma_disp", "#beta^{-} clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiBeta_ThreeSigma", "#beta^{+} clusters", nPt,ptMin,ptMax,1000,0.,100.));
  fOutputContainer->Add(new TH2F("hTracksOfAntiBeta_ThreeSigma_disp", "#beta^{+} clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));

  fOutputContainer->Add(new TH1F("hTracksOfOthers", "Other clusters", nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hTracksOfOthers_disp", "Other clusters, Disp cut", nPt,ptMin,ptMax,1000,0.,100.));
  

//=================

  fOutputContainer2->Add(new TH2F("hTracksOfPi_OneSigma_label", "PDG #pi^{+}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfPi_OneSigma_disp_label", "PDG #pi^{+}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPi_OneSigma_label", "PDG #pi^{-}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPi_OneSigma_disp_label", "PDG #pi^{-}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));

  fOutputContainer2->Add(new TH2F("hTracksOfPi_ThreeSigma_label", "PDG #pi^{+}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfPi_ThreeSigma_disp_label", "PDG #pi^{+}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPi_ThreeSigma_label", "PDG #pi^{-}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPi_ThreeSigma_disp_label", "PDG #pi^{-}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));


  fOutputContainer2->Add(new TH2F("hTracksOfPr_OneSigma_label", "PDG p-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfPr_OneSigma_disp_label", "PDG p-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPr_OneSigma_label", "PDG #bar{p}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPr_OneSigma_disp_label", "PDG #bar{p}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));

  fOutputContainer2->Add(new TH2F("hTracksOfPr_ThreeSigma_label", "PDG p-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfPr_ThreeSigma_disp_label", "PDG p-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPr_ThreeSigma_label", "PDG #bar{p}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiPr_ThreeSigma_disp_label", "PDG #bar{p}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));


  fOutputContainer2->Add(new TH2F("hTracksOfKa_OneSigma_label", "PDG K^{+}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfKa_OneSigma_disp_label", "PDG K^{+}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiKa_OneSigma_label", "PDG K^{-}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiKa_OneSigma_disp_label", "PDG K^{-}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));

  fOutputContainer2->Add(new TH2F("hTracksOfKa_ThreeSigma_label", "PDG K^{+}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfKa_ThreeSigma_disp_label", "PDG Ka^{+}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiKa_ThreeSigma_label", "PDG K^{-}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiKa_ThreeSigma_disp_label", "PDG K^{-}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));


  fOutputContainer2->Add(new TH2F("hTracksOfBeta_OneSigma_label", "PDG #beta^{-}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfBeta_OneSigma_disp_label", "PDG #beta^{-}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiBeta_OneSigma_label", "PDG #beta^{+}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiBeta_OneSigma_disp_label", "PDG #beta^{+}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));

  fOutputContainer2->Add(new TH2F("hTracksOfBeta_ThreeSigma_label", "PDG #beta^{-}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfBeta_ThreeSigma_disp_label", "PDG #beta^{-}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiBeta_ThreeSigma_label", "PDG #beta^{+}-s vs p_{T}", nPt,ptMin,ptMax, 8000, -4000, 4000));
  fOutputContainer2->Add(new TH2F("hTracksOfAntiBeta_ThreeSigma_disp_label", "PDG #beta^{+}-s vs p_{T}, dispersion cut", nPt,ptMin,ptMax, 8000, -4000, 400));


//=================


    fOutputContainer2->Add(new TH1F("hTracksOfPiClose","Tracks of #pi^{#pm} close to cluster", nPt,ptMin,ptMax));
    fOutputContainer2->Add(new TH1F("hTracksOfPiCloseDispOK","Tracks of #pi^{#pm} close to cluster", nPt,ptMin,ptMax));
    fOutputContainer2->Add(new TH1F("hTracksOfPrClose","Tracks of p#bar{p} close to cluster", nPt,ptMin,ptMax));
    fOutputContainer2->Add(new TH1F("hTracksOfPrCloseDispOK","Tracks of p#bar{p} close to cluster", nPt,ptMin,ptMax));
    fOutputContainer2->Add(new TH1F("hTracksOfKaClose","Tracks of K^{#pm} close to cluster", nPt,ptMin,ptMax));
    fOutputContainer2->Add(new TH1F("hTracksOfKaCloseDispOK","Tracks of Ka^{#pm} close to cluster", nPt,ptMin,ptMax));
    fOutputContainer2->Add(new TH1F("hTracksOfOthersClose","Tracks of K^{#pm} close to cluster", nPt,ptMin,ptMax));
    fOutputContainer2->Add(new TH1F("hTracksOfOthersCloseDispOK","Tracks of Ka^{#pm} close to cluster", nPt,ptMin,ptMax));
        
    
   fOutputContainer2->Add(new TH1F("hTracks_matched","Nr of matched tracks per cluster",20,0,20));

   fOutputContainer2->Add(new TH2F("hEnTrackvsClust","Energies of track vs cluster",nPt,ptMin,ptMax,nPt,ptMin,ptMax));

   fOutputContainer2->Add(new TH1F("hEmcCPVDistance","EMS to CPV distance",100,0.,100.));

   fOutputContainer2->Add(new TH1F("hDistance","Distance from cluster to associated track, cm",100,0.,100.));

// fOutputContainer2->Add(new TH2F("hTracks_matched_PID_NC",
  
//======================== MC histograms!!!!  ============================================

  fOutputContainer2->Add(new TH2F("hMCPdgvsPt","Pdg code vs particle p_{T}",nPt,ptMin,ptMax, 8000,-4000,4000));
  fOutputContainer2->Add(new TH1F("hMCPt","MC particles p_{T}",nPt,ptMin,ptMax));


  fOutputContainer2->Add(new TH2F("fhPi0MC", "MC distribution of #pi^{0}-s ", nPt,ptMin,ptMax, 240, -1.2, 1.2));
  Sumw2Histogram("fhPi0MC");
  fOutputContainer2->Add(new  TH2F("fhEtaMC", "MC distribution of #eta^{0}-s ", nPt,ptMin,ptMax, 240, -1.2, 1.2));
  Sumw2Histogram("fhEtaMC");
  fOutputContainer2->Add(new  TH2F("fhEtaPrimeMC", "MC distribution of #eta'", nPt,ptMin,ptMax,240, -1.2, 1.2));
  Sumw2Histogram("fhEtaPrimeMC");
  fOutputContainer2->Add(new  TH2F("fhOmegaMC", "MC distribution of #omega", nPt,ptMin,ptMax, 240, -1.2, 1.2));
  Sumw2Histogram("fhOmegaMC");
  fOutputContainer2->Add(new  TH2F("fhK0SMC", "MC distribution of #K_{0}^{S}-s", nPt,ptMin,ptMax, 240, -1.2, 1.2));
    Sumw2Histogram("fhK0SMC");
  fOutputContainer2->Add(new  TH2F("fhK0LMC", "MC distribution of #K_{0}^{L}-s", nPt,ptMin,ptMax, 240, -1.2, 1.2));
  Sumw2Histogram("fhK0LMC");
  fOutputContainer2->Add(new  TH2F("fhGammaMC_all", "MC distribution of #gamma-s", nPt,ptMin,ptMax, 240, -1.2, 1.2));
  Sumw2Histogram("fhGammaMC_all");
  fOutputContainer2->Add(new  TH2F("fhGammaMC_true", "MC distribution of #gamma-s", nPt,ptMin,ptMax, 240, -1.2, 1.2));
  Sumw2Histogram("fhGammaMC_true");
  fOutputContainer2->Add(new  TH2F("fhGammaMCSources", "MC distribution of #gamma-s", nPt,ptMin,ptMax, 8000, -4000,4000));
  Sumw2Histogram("fhGammaMCSources");
    fOutputContainer2->Add(new  TH2F("fhProtonMC", "MC distribution of protons", nPt,ptMin,ptMax, 240, -1.2, 1.2));
  Sumw2Histogram("fhProtonMC");
  fOutputContainer2->Add(new  TH2F("fhNeutrMC", "MC distribution of neutrons", nPt,ptMin,ptMax, 240, -1.2, 1.2));
  Sumw2Histogram("fhNeutrMC");
  fOutputContainer2->Add(new  TH2F("fhKchMC", "MC distribution of neutrons", nPt,ptMin,ptMax, 240, -1.2, 1.2));
  Sumw2Histogram("fhKchMC");
  fOutputContainer2->Add(new  TH2F("fhPdgvsPt_MCTracks","MC Particle PDG code vs Pt", nPt,ptMin,ptMax , 8000, -4000, 4000));
  Sumw2Histogram("fhPdgvsPt_MCTracks");
  fOutputContainer2->Add(new  TH2F("fhPrimMC ","MC Particle PDG code vs Pt, primary particles only", nPt,ptMin,ptMax , 40, -20, 20));

  fOutputContainer2->Add(new TH2F("fhGenvsMeas","Generated vs measured energy",nPt,ptMin,ptMax,nPt,ptMin,ptMax ));
  fOutputContainer2->Add(new TH2F("fhGenvsMeas_corr","Generated vs measured energy",nPt,ptMin,ptMax,nPt,ptMin,ptMax ));

   for(Int_t iCut = 0; iCut < Ncut; iCut++)
   {
    fOutputContainer2->Add(new TH2F(Form("hClustPdgvsPt_%s", cut[iCut].Data()),"Cluster pdg vs p_{T}." ,nPt,ptMin,ptMax,8000, -4000,4000));
    fOutputContainer2->Add(new TH2F(Form("hClustPdgvsPt_%s_naive", cut[iCut].Data()),"Cluster pdg vs p_{T} (naive)." ,nPt,ptMin,ptMax,8000, -4000,4000));
    Sumw2Histogram(Form("hClustPdgvsPt_%s", cut[iCut].Data()));
    Sumw2Histogram(Form("hClustPdgvsPt_%s_naive", cut[iCut].Data()));
    fOutputContainer2->Add(new TH2F(Form("hMatrixEff_%s", cut[iCut].Data()), "Efficiency matrix ", 400, 0., 40., 400, 0., 40.)); 
   }

   fOutputContainer->Add(new TH1F("htest","Count events", 1, 0, 1));
   fOutputContainer->Add(new TH1F("hEventCounter","Count accepted events", 1, 0, 1));
   fOutputContainer->Add(new TH1F("hEventCounterCentrality","Count events in centrality bins", 8, 0. , 8.));
   fOutputContainer->Add(new TH1F("hPHOSCounter","Count accepted photons", 2, 0, 2));
   fOutputContainer->Add(new TH1F("hPHOSvsEMCAL", "PHOS and EMCAL clusters", 2, 0, 2));
   fOutputContainer->Add(new TH1I("hPhosChargedNeutral", "PHOS clusters 1-all, 2-neutral, 3-charged", 3, 0, 3));
   fOutputContainer2->Add(new TH1F("htestmy","htestmy",400,0,40));
   fOutputContainer2->Add(new TH1F("hteststandard","hteststandard",400,0,40));
   fOutputContainer2->Add(new TH1I("hpid3", "Pid #sigma<3, #pi, p, K, #beta, Undef", 5, 0., 5.));
   fOutputContainer2->Add(new TH1I("hpid1", "Pid #sigma<1, #pi, p, K, #beta, Undef", 5, 0., 5.));
   fOutputContainer2->Add(new TH1F("hEventCounterMC","Count events",1,0,1));

  PostData(1, fOutputContainer);
  PostData(2, fOutputContainer2);
}

/*
================================================================================
  // Main loop, called for each event

================================================================================
*/
void AliAnalysisTaskGammaPHOSPP::UserExec(Option_t *) 
{

  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) 
  {
     Printf("ERROR: Could not retrieve event");
     return;
  }

  fPHOSGeo = AliPHOSGeometry::GetInstance();

  if (fPHOSGeo)
  {
      AliInfo("PHOS geometry not initialized, initializing it for you");

      if(fEvent->GetRunNumber() < 224994)
        fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP"); // Run1 geometry
      else
        fPHOSGeo = AliPHOSGeometry::GetInstance("Run2");
  }


  if(fEvent->GetRunNumber() > 224994)
    fWeightFunction= new TF1("fWeightFunction", "1.0", 0., 99999.) ;
  else
    fWeightFunction= new TF1("fWeightFunction", "(3.02640e-01 + 8.59672e-01*x + 5.39777e-01*x*x)/(1.+ 9.82832e-02 *x+1.27487e+00 *x*x)+3.73416e-03*x", 0., 99999.) ;

  fAllEventCounter++;

/*----------------------------Vertex------------------------------------------*/

  const AliAODVertex *esdVertex5 =    fEvent->GetPrimaryVertex();
  const AliAODVertex *esdVertexSPD  = fEvent->GetPrimaryVertexSPD();
  
  Double_t fVtx0[3] = {0, 0, 0};
  Double_t fVtx5[3] = {esdVertex5->GetX(), esdVertex5->GetY(), esdVertex5->GetZ()};


/*--------------------------Filter events-------------------------------------*/

  Bool_t acceptEvent = AcceptEvent(fEvent);
  if(!acceptEvent)
  {
      Printf("Event nr %i is rejected", fAllEventCounter);
      return;  
  }
  else
      Printf("Event nr %i is accepted", fAllEventCounter);

/*--------------------------PHOS event----------------------------------------*/

  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton", 50) ;
  
/*-------------------------PID------------------------------------------------*/
    
fPIDResponse -> SetUseTPCMultiplicityCorrection(kFALSE);
AliPIDCombined *pidcomb=new AliPIDCombined();
pidcomb->SetDefaultTPCPriors();
//pidcomb->SetEnablePriors(kFALSE);
pidcomb->SetSelectedSpecies(AliPID::kSPECIESC);
pidcomb->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF|AliPIDResponse::kDetITS|AliPIDResponse::kDetTRD);


/*----------------------------------------------------------------------------*/
/* Event centrality */
  //always zero centrality
  fEventCentrality = 0;
  
  fEventCentrality = GetEventCentrality(fEvent);

/*----------------------------------------------------------------------------*/
/* Process MC */
      
  fMCArray = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());
  
  Notify();

  ProcessMC() ;

/*----------------------------------------------------------------------------*/
/* Count PHOS and EMCAL clusters */

  PHOSvsEMCALClusters();
  
/*----------------------------------------------------------------------------*/
/* PHOS cells and clusters */

   AnalyzeCells();

   fInPHOS = 0 ;
   
   for(Int_t ic = 0; ic < fEvent->GetNumberOfCaloClusters(); ic++)
   {
      AliAODCaloCluster *clu1 = fEvent->GetCaloCluster(ic);
      SelectClusters(clu1);
   }
   
   for(Int_t iph = 0; iph < fInPHOS; iph++)
   {
     AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(iph) ;
     FillOnePhotonHistograms(ph);
   }

   FillTwoPhotonHistograms();
  
   MixPhotons();

/*----------------------------------------------------------------------------*/

  FillHistogram("hEventCounter",0.5);
  FillHistogram("hEventCounterCentrality", fEventCentrality + 0.5);
  fEventCounter++;
  
/*----------------------------------------------------------------------------*/  
// Post output data.
 
  PostData(1, fOutputContainer);
  PostData(2, fOutputContainer2);
}

/*----------------------------------------------------------------------------*/
/*============================================================================*/
/*----------------------------------------------------------------------------*/

void AliAnalysisTaskGammaPHOSPP::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
/*

  Int_t nPP = fnCINT1B - fnCINT1A - fnCINT1C + 2*fnCINT1E;
  FillHistogram("hTrigClass",1,fnCINT1B);
  FillHistogram("hTrigClass",2,fnCINT1A);
  FillHistogram("hTrigClass",3,fnCINT1C);
  FillHistogram("hTrigClass",4,fnCINT1E);
  FillHistogram("hTrigClass",5,nPP);
  Printf("fnCINT1B=%d, fnCINT1A=%d ,fnCINT1C=%d ,fnCINT1E=%d, nPP=%d",
     fnCINT1B,fnCINT1A,fnCINT1C,fnCINT1E,nPP);
  */ 
}

//________________________________________________________________________

Bool_t AliAnalysisTaskGammaPHOSPP::AcceptEvent(AliAODEvent *aodEvent)
{

  FillHistogram("hSelEvents",0) ; // All events accepted by Physics Selection

  TString trigClasses = aodEvent->GetFiredTriggerClasses();
  Printf("Event nr %i with triggers %s, period %i.", fAllEventCounter, 
                                                     trigClasses.Data(), 
                                                     fEvent->GetPeriodNumber());
                                                     
  if (trigClasses.Contains("FAST")  && !trigClasses.Contains("ALL")) 
  {
    AliWarning(Form("Skip event with triggers %s",trigClasses.Data()));
    return kFALSE;
  }

  if (trigClasses.Contains("CINT1B")) fnCINT1B++;
  if (trigClasses.Contains("CINT1A")) fnCINT1A++;
  if (trigClasses.Contains("CINT1C")) fnCINT1C++;
  if (trigClasses.Contains("CINT1-E")) fnCINT1E++;

  // Event selection flags

  Bool_t eventVtxExist    = kFALSE;
  Bool_t eventVtxZ10cm    = kFALSE;
  Bool_t eventPileup      = kFALSE;
  Bool_t eventV0AND       = kFALSE;


  Int_t eventNumberInFile = aodEvent->GetEventNumberInFile();

  if (aodEvent->GetPrimaryVertex()->GetNContributors() <1 && !fMCArray)
      eventVtxExist    = kFALSE;
    else     
      eventVtxExist    = kTRUE; 
   

  if (aodEvent->IsPileupFromSPD())
    eventPileup = kTRUE;

  const AliAODVertex *esdVertex5 =    aodEvent->GetPrimaryVertex();
  const AliAODVertex *esdVertexSPD  = aodEvent->GetPrimaryVertexSPD();

 FillHistogram("hNvertexTracks", esdVertex5->GetNContributors());


  if (aodEvent->GetPrimaryVertex() && aodEvent->GetPrimaryVertex()->GetNContributors() > 0)
  {
     FillHistogram("hZvertex", esdVertex5->GetZ());
     if (TMath::Abs(esdVertex5->GetZ()) < 10.)
	eventVtxZ10cm = kTRUE;
  }

  if (aodEvent->IsPileupFromSPD()) 
    {
       eventPileup = kTRUE;
         /*
         TClonesArray *pileupVertices = event->GetPileupVerticesSPD();
         Int_t nPileupVertices = pileupVertices->GetEntriesFast();
         */
       FillHistogram("hNPileupVtx", aodEvent->GetNumberOfPileupVerticesSPD());
         for (Int_t puVtx = 0;  puVtx < fEvent->GetNumberOfPileupVerticesSPD(); puVtx++) 
         {
           Double_t dZpileup = esdVertexSPD->GetZ() - fEvent->GetPileupVertexSPD(puVtx)->GetZ();
           FillHistogram("hZPileupVtx",dZpileup);
         }
    }

  ULong64_t trigmask = aodEvent->GetTriggerMask();
  
  eventV0AND = (trigmask & (1ull << (AliTriggerAnalysis::kV0AND-1)));

  // Fill event statistics for different selection criteria

    FillHistogram("hSelEvents",1) ;
  if (eventVtxExist) 
    FillHistogram("hSelEvents",2) ;
  if (eventVtxExist && eventVtxZ10cm)
    FillHistogram("hSelEvents",3) ;
  if (eventVtxExist && eventVtxZ10cm && eventV0AND)
    FillHistogram("hSelEvents",4) ;
  if (eventVtxExist && eventVtxZ10cm && eventV0AND && eventPileup)
    FillHistogram("hSelEvents",5) ;
  if (eventPileup)
    FillHistogram("hSelEvents",6) ;
  if(eventV0AND)
    FillHistogram("hSelEvents",7) ;
  if(eventVtxZ10cm)
    FillHistogram("hSelEvents",8) ;  
        
  if(esdVertex5->GetNContributors() < 1  && !fMCArray)
     return kFALSE;
  if(esdVertexSPD->GetNContributors() < 1 && !fMCArray) 
     return kFALSE;
  if (!eventVtxExist) return kFALSE;
  if (!eventVtxZ10cm) return kFALSE;
  if (eventPileup)    return kFALSE;
  
  return kTRUE;

}

/*============================================================================*/
Int_t AliAnalysisTaskGammaPHOSPP::GetEventCentrality(AliAODEvent *event)
{


  //Calculate charged multiplicity   Bool_t   IsCPVOK(void)const {return fCpv;}
  
  for (Int_t i=0; i < event->GetNumberOfTracks(); i++) 
  {
      AliAODTrack* track = dynamic_cast<AliAODTrack*>(event->GetTrack(i)) ;
      AliVTrack*   ttrack=(AliVTrack*)event->GetTrack(i);

      FillHistogram("hTrackCharge", track->Charge());
  }

  Int_t trackMult = event->GetNumberOfTracks() ;

  Float_t tV0A = event->GetVZEROData()->GetV0ATime();
  Float_t tV0C = event->GetVZEROData()->GetV0CTime();
  FillHistogram("hV0Atime", tV0A);
  FillHistogram("hV0Atime", tV0C);
  FillHistogram("hV0AV0Ctime", tV0A,tV0C);  
  FillHistogram("hTrackMult", trackMult+0.5) ;

  Int_t centr = 0;
  if(trackMult <= 2)
    centr = 0 ;
  else 
    if(trackMult <= 5)
      centr=1 ;
    else
      if(trackMult<=9)
        centr=2 ;
      else
        if(trackMult<=14)
          centr=3 ;
        else
          if(trackMult<=22)
            centr=4 ;
          else
            if(trackMult<=35)
              centr=5 ;
            else
              if(trackMult<=50)
                centr=6 ;
              else
                centr=7 ;

  return centr;

}

//_____________________________________________________________________________

void AliAnalysisTaskGammaPHOSPP::FillHistogram(const char * key,Double_t x)const
{
  //FillHistogram
  TH1I * tmpI = dynamic_cast<TH1I*>((fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key)) ) ;

  if(tmpI)
  {
    tmpI->Fill(x) ;
    return ;
  }
  TH1F * tmpF = dynamic_cast<TH1F*>((fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key)) ) ;

  if(tmpF)
  {
    tmpF->Fill(x) ;
    return ;
  }
  TH1D * tmpD = dynamic_cast<TH1D*>((fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key)) ) ;

  if(tmpD)
  {
    tmpD->Fill(x) ;
    return ;
  }
  AliInfo(Form("can not find histogram <%s> (Fill I)",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::FillHistogram(const char * key,Double_t x,Double_t y)const
{
  //FillHistogram
  TObject * tmp = (fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key)) ;
  if(!tmp)
  {
    AliInfo(Form("can not find histogram <%s> (Fill II) ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1F"))
  {
    ((TH1F*)tmp)->Fill(x,y) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F"))
  {
    ((TH2F*)tmp)->Fill(x,y) ;
    return ;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName())) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const
{
  //Fills 1D histograms with key
  TObject * tmp = (fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key)) ;

  if(!tmp)

  {
    AliInfo(Form("can not find histogram <%s> (Fill III) ",key)) ;
    return ;
  }

  if(tmp->IsA() == TClass::GetClass("TH2F"))
  {
    ((TH2F*)tmp)->Fill(x,y,z) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F"))
  {
    ((TH3F*)tmp)->Fill(x,y,z) ;
    return ;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPP::Sumw2Histogram(const char * key) const
{
  //Fills 1D histograms with key
  TObject * tmp = (fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key)) ;
  if(!tmp)
  {
    AliInfo(Form("can not find histogram <%s> (Sumw2) ",key)) ;
    return ;
  }
    if(tmp->IsA() == TClass::GetClass("TH1F"))
  {
    ((TH1F*)tmp)->Sumw2() ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F"))
  {
    ((TH2F*)tmp)->Sumw2() ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F"))
  {
    ((TH3F*)tmp)->Sumw2() ;
    return ;
  }
}

//--------------------------------------------------------------------------------
void AliAnalysisTaskGammaPHOSPP::ProcessMC()
{
  if(!fMCArray) 
   return; 
  char partName[10] ;

  FillHistogram("hEventCounterMC",0.5);
   
   for(Int_t i = 0; i < fMCArray->GetEntriesFast(); i++)
   {
       AliAODMCParticle* particle =  (AliAODMCParticle*) fMCArray->At(i);
       FillHistogram("hMCTrackCharge", particle->Charge());

       if(TMath::Abs(particle->GetPdgCode()) == 111)
       {
          FillHistogram("fhPi0MC", particle->Pt(), particle->Y(), Weight(particle));
          snprintf(partName, 10, "pi0") ;
       }
       else if(TMath::Abs(particle->GetPdgCode()) == 221)
           FillHistogram("fhEtaMC", particle->Pt(), particle->Y(), Weight(particle)); 
       else if(TMath::Abs(particle->GetPdgCode()) == 331)
           FillHistogram("fhEtaPrimeMC", particle->Pt(), particle->Y(), Weight(particle)); 
       else if(TMath::Abs(particle->GetPdgCode()) == 223)
           FillHistogram("fhOmegaMC", particle->Pt(), particle->Y(), Weight(particle)); 
       else if(TMath::Abs(particle->GetPdgCode()) == 130)
           FillHistogram("fhK0LMC", particle->Pt(), particle->Y(), Weight(particle));        
       else if(TMath::Abs(particle->GetPdgCode()) == 310)
           FillHistogram("fhK0SMC", particle->Pt(), particle->Y(),  Weight(particle));         
       else if(TMath::Abs(particle->GetPdgCode() == 2212)) 
           FillHistogram("fhProtonMC",particle->Pt(), particle->Y(), Weight(particle));
       else if(TMath::Abs(particle->GetPdgCode()==2112)) 
           FillHistogram("fhNeutrMC",particle->Pt(), particle->Y(), Weight(particle));
       else if(TMath::Abs(particle->GetPdgCode()==321)) 
           FillHistogram("fhKchMC",particle->Pt(), particle->Y(), Weight(particle));
       else
       if(TMath::Abs(particle->GetPdgCode()) == 22)
       {
            FillHistogram("fhGammaMC_all", particle->Pt(), particle->Y(), Weight(particle)); 
            Int_t label1 = particle->GetMother(); 
            if(label1==-1) 
            { 	     
 	       FillHistogram("fhGammaMC_true", particle->Pt(), particle->Y(), Weight(particle));
	       FillHistogram("fhGammaMCSources", particle->Pt(), 22, Weight(particle));
             }
             else
             {
 		    AliAODMCParticle *particle1 =  (AliAODMCParticle*) fMCArray->At(label1);
	            if(particle1->GetPdgCode()!=111)
                    {		
                        FillHistogram("fhGammaMCSources", particle->Pt(), 
                        particle1->GetPdgCode(), Weight(particle));
			FillHistogram("fhGammaMC_true",particle->Pt(),particle->Y(), Weight(particle));
                     }
	             else
                     {
			 if(particle1->IsPrimary())
                         {
			    FillHistogram("fhGammaMC_true",particle->Pt(), particle->Y(), Weight(particle));
                            FillHistogram("fhGammaMCSources",particle->Pt(),111,Weight(particle));
                          }
			  else
                          {
			      Int_t label2=particle1->GetMother();
			      AliAODMCParticle *particle2 =  (AliAODMCParticle*) fMCArray->At(label2);
			      if(TMath::Abs(particle2->GetPdgCode())!=130 
			                  && TMath::Abs(particle2->GetPdgCode())!=310)
                              {
			         FillHistogram("fhGammaMC_true",particle->Pt(),
                                                particle->Y(), Weight(particle));
			         FillHistogram("fhGammaMCSources", particle->Pt(),111);
		              }
			      else 
			         FillHistogram("fhGammaMCSources", particle->Pt(),
				                  particle2->GetPdgCode(), Weight(particle));
			  }  
		     }
              }
       }      
       continue ;
   }
  
  
}
//------------------------------------------------------------------------------
void AliAnalysisTaskGammaPHOSPP::AnalyzeCells()
{
  AliAODCaloCells *cells = fEvent->GetPHOSCells();

  Int_t multCells = cells->GetNumberOfCells();
  FillHistogram("hCellMultEvent",multCells);

  Float_t  energy;
  Int_t    mod1, relId[4], cellAbsId, cellX, cellZ;

  // Single loop over cells

  Int_t nCellModule[4] = {0, 0, 0, 0};
  for (Int_t iCell=0; iCell<multCells; iCell++) 
  {
    cellAbsId = cells->GetCellNumber(iCell);
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    mod1  = relId[0];
    cellX = relId[2];
    cellZ = relId[3] ;
    energy = cells->GetAmplitude(iCell);
    FillHistogram("hCellEnergy",energy);
    if      (mod1==1) 
    {
      nCellModule[0]++;
      FillHistogram("hCellEnergyM1",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM1",cellX,cellZ,1.);
      FillHistogram("hCellEXZM1",cellX,cellZ,energy);
     }
    else if (mod1==2) 
    {
      nCellModule[1]++;
      FillHistogram("hCellEnergyM2",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM2",cellX,cellZ,1.);
      FillHistogram("hCellEXZM2",cellX,cellZ,energy);
    }
    else if (mod1==3) 
    {
      nCellModule[2]++;
      FillHistogram("hCellEnergyM3",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM3",cellX,cellZ,1.);
      FillHistogram("hCellEXZM3",cellX,cellZ,energy);
    }
    else if (mod1==4) 
    {
      nCellModule[3]++;
      FillHistogram("hCellEnergyM4",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM4",cellX,cellZ,1.);
      FillHistogram("hCellEXZM4",cellX,cellZ,energy);
    }
  }
  FillHistogram("hCellMultEventM1",nCellModule[0]);
  FillHistogram("hCellMultEventM2",nCellModule[1]);
  FillHistogram("hCellMultEventM3",nCellModule[2]);
  FillHistogram("hCellMultEventM4",nCellModule[3]);
}

//------------------------------------------------------------------------------

void AliAnalysisTaskGammaPHOSPP::SelectClusters(AliAODCaloCluster *clu1)
{
  //Analyze clusters and select photons for analysis

  Int_t  multPHOSClust[5]  = {0, 0, 0, 0, 0}; 
  
  TLorentzVector p1, p11;
  
  Int_t    mod1, relId[4], cellAbsId, cellX, cellZ;
  Float_t  position[3];  
  Int_t digMult;
  Double_t energy;
  Double_t weight = 1.0;



//  for (Int_t ic = 0; ic < multClust; ic ++) 
//  {
//    AliAODCaloCluster *clu1 = fEvent->GetCaloCluster(ic);

    if(!clu1->IsPHOS() ) return;
         
    clu1->GetPosition(position);
    TVector3 global1(position) ;
    //fPHOSGeo->GlobalPos2RelId(global1,relId) ;
    
    cellAbsId = clu1->GetCellAbsId(0);
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    mod1  = relId[0] ;
    cellX = relId[2] ;
    cellZ = relId[3] ;
    energy = clu1->E();
    digMult = clu1->GetNCells(); 
      
    if (mod1 < 1 || mod1 > 4) 
    {
      Printf("Wrong module number %d", mod1);
      return;
    }
    
    if(clu1->GetType() == AliVCluster::kPHOSCharged)  
      return;
    if(fEvent->GetRunNumber() > 224994 && !fMCArray && TMath::Abs(clu1->GetTOF()) > 12.5e-9) 
      return; // TOF cut for real data only!
    
    multPHOSClust[0]++;
    FillHistogram("hClusterEnergy",energy);
    FillHistogram("hCellMultClu_all",digMult);
    FillHistogram("hClusterEvsN_all",energy,digMult);
    FillHistogram(Form("hClusterEnergyM%d", mod1), energy);
    FillHistogram(Form("hClusterEvsN_all_M%d", mod1), energy, digMult);
    FillHistogram(Form("hCellMultClu_all_M%d", mod1), digMult);
    FillHistogram(Form("hCluNXZM%d", mod1), cellX, cellZ, 1.);
    if(clu1->GetEmcCpvDistance() > 2.5)
    {
        FillHistogram("hClusterEvsN_cpv",energy,digMult);
        FillHistogram("hCellMultClu_cpv",digMult);
        FillHistogram(Form("hClusterEvsN_cpv_M%d", mod1), energy, digMult);
        FillHistogram(Form("hCellMultClu_cpv_M%d", mod1), digMult);
    }
    if(clu1->Chi2() < 2.5)
    {
        FillHistogram("hClusterEvsN_disp",energy,digMult);
        FillHistogram("hCellMultClu_disp",digMult);
        FillHistogram(Form("hClusterEvsN_disp_M%d", mod1), energy, digMult);
        FillHistogram(Form("hCellMultClu_disp_M%d", mod1), digMult);
    }
    if( clu1->GetEmcCpvDistance() > 2.5 && clu1->Chi2() < 2.5)
    {
        FillHistogram("hClusterEvsN_both",energy,digMult);
        FillHistogram("hCellMultClu_both",digMult);
        FillHistogram(Form("hClusterEvsN_both_M%d", mod1), energy,digMult);
        FillHistogram(Form("hCellMultClu_both_M%d", mod1), digMult);
    }

    if (mod1==1) 
      multPHOSClust[1]++;
    else if (mod1==2) 
      multPHOSClust[2]++;
    else if (mod1==3) 
      multPHOSClust[3]++;
    else if (mod1==4) 
      multPHOSClust[4]++;    
     
     AliPHOSAodCluster cluPHOS1(*clu1);

     cluPHOS1.GetMomentum(p1 , fVtx0);
     cluPHOS1.GetMomentum(p11, fVtx5);
     
      Double_t pAbs = p1.P();
      Double_t pT   = p1.Pt();
      Double_t pX   = p1.Px();
      Double_t pY   = p1.Py();
      if (pAbs<1.e-10) return;
      Double_t kappa = pAbs - TMath::Power(0.135, 2)/4./pAbs;
   
      FillHistogram("hPhotonKappa", kappa);
      FillHistogram("hPhotonPt", pT);
      FillHistogram("hPhotonPx", pX);
      FillHistogram("hPhotonPy", pY);

      if(clu1->E()       < 0.3) return;
      if(clu1->GetNCells() < 3) return ;       
      if(clu1->GetM02() < 0.2)  return ;    

      FillHistogram("hEmcCPVDistance", clu1->GetEmcCpvDistance());
      TestMatchingTrackPID(clu1, p11.Pt());
     
      new((*fPHOSEvent)[fInPHOS]) AliCaloPhoton(p1.X(),p1.Py(),p1.Z(),p1.E()) ;
      AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(fInPHOS) ;
      ph->SetModule(mod1) ;
      ph->SetMomV2(&p11) ;
      ph->SetNCells(clu1->GetNCells());
      ph->SetEMCx(global1.X()); 
      ph->SetEMCy(global1.Y());
      ph->SetEMCz(global1.Z());
      ph->SetDispBit(clu1->Chi2() < 2.5*2.5);
      ph->SetCPVBit(clu1->GetEmcCpvDistance() > 2.5);
      ph->SetBC(TestBC(clu1->GetTOF()));
      ph->SetPrimary(GetPrimaryLabel(clu1));
      ph->SetPrimaryAtVertex(GetPrimaryLabelAtVertex(clu1));
      ph->SetWeight(weight); 
      
      FillHistogram("hvt0vsvt5", p11.Pt()- p1.Pt());
      FillHistogram("hBC", TestBC(clu1->GetTOF()) + 0.5);
      FillHistogram("hTOF", clu1->GetTOF());
      FillHistogram("hWeights", ph->GetWeight());      

      fInPHOS++ ;
      
//   }
/*     
   FillHistogram("hPHOSClusterMult",   multPHOSClust[0]);
   FillHistogram("hPHOSClusterMultM1", multPHOSClust[1]);
   FillHistogram("hPHOSClusterMultM2", multPHOSClust[2]);
   FillHistogram("hPHOSClusterMultM3", multPHOSClust[3]);
   FillHistogram("hPHOSClusterMultM4", multPHOSClust[4]);
*/
}

//===========================================================================

void AliAnalysisTaskGammaPHOSPP::FillOnePhotonHistograms(AliCaloPhoton *ph)
{

      TLorentzVector p11 = *(ph->GetMomV2());
      Double_t pT = p11.Pt();
      Double_t weight = ph->GetWeight();
      Int_t pdg = 0;
      Int_t pdg_naive = 0;
      if(fMCArray)
      {
        pdg = ((AliAODMCParticle*)fMCArray->At(ph->GetPrimaryAtVertex())) -> GetPdgCode();
        pdg_naive = ((AliAODMCParticle*)fMCArray->At(ph->GetPrimary())) -> GetPdgCode();
        weight = Weight((AliAODMCParticle*)fMCArray->At(ph->GetPrimaryAtVertex()));
        Printf("pdg_naive = %d, pdg = %d", pdg_naive, pdg);
      }
      Int_t mod1 = ph->Module();

      Double_t ww = 1;
      if(pdg == 22) FillHistogram("hMatrixEff_all", p11.Pt(), TestGammaPt(ph), weight);       

      FillHistogram("hClustPt_all", p11.Pt(), weight );
      FillHistogram("hClustPdgvsPt_all", p11.Pt(), pdg, weight);
      FillHistogram("hClustPdgvsPt_all_naive",  p11.Pt(), pdg_naive, weight);

      FillHistogram("hClustPtvsNcl_all", p11.Pt(), ph->GetNCells() + 0.5, weight);
      FillHistogram("hCentralityvsClustPt_all", p11.Pt(), fEventCentrality + 0.5, weight);	
      FillHistogram(Form("hClustPt_all_mod%d", mod1),  p11.Pt(), weight);

      if(ph->IsCPVOK())
      {
        FillHistogram("hClustPt_cpv", p11.Pt(), weight );
        FillHistogram("hClustPdgvsPt_cpv", p11.Pt(), pdg, weight); 
        FillHistogram("hClustPdgvsPt_cpv_naive",  p11.Pt(), pdg_naive, weight);  

        FillHistogram("hClustPtvsNcl_cpv", p11.Pt(), ph->GetNCells()  + 0.5, weight);
        FillHistogram("hCentralityvsClustPt_cpv", p11.Pt(), fEventCentrality + 0.5, weight);	
        FillHistogram(Form("hClustPt_cpv_mod%d", mod1),  p11.Pt(), weight);

        if(pdg == 22) FillHistogram("hMatrixEff_cpv", p11.Pt(), TestGammaPt(ph), weight);       

      }
      if(ph->IsDispOK())
      {
        FillHistogram("hClustPt_disp", p11.Pt(), weight );
        FillHistogram("hClustPdgvsPt_disp", p11.Pt(), pdg, weight);
        FillHistogram("hClustPdgvsPt_disp_naive",  p11.Pt(), pdg_naive, weight);

        FillHistogram("hClustPtvsNcl_disp", p11.Pt(), ph->GetNCells()  + 0.5, weight);
        FillHistogram("hCentralityvsClustPt_disp", p11.Pt(), fEventCentrality + 0.5, weight);	
        FillHistogram(Form("hClustPt_disp_mod%d", mod1),  p11.Pt(), weight);

        if(pdg == 22) FillHistogram("hMatrixEff_disp", p11.Pt(), TestGammaPt(ph), weight);       

      }
      if( ph->IsCPVOK() && ph->IsDispOK())
      { 
        FillHistogram("hClustPt_both", p11.Pt(), weight );
        FillHistogram("hClustPdgvsPt_both", p11.Pt(), pdg, weight);
        FillHistogram("hClustPdgvsPt_both_naive",  p11.Pt(), pdg_naive, weight);

        FillHistogram("hClustPtvsNcl_both", p11.Pt(), ph->GetNCells() + 0.5, weight);
        FillHistogram("hCentralityvsClustPt_both", p11.Pt(), fEventCentrality + 0.5, weight);	
        FillHistogram(Form("hClustPt_both_mod%d", mod1),  p11.Pt(), weight);

        if(pdg == 22) FillHistogram("hMatrixEff_both", p11.Pt(), TestGammaPt(ph), weight);       
      }
}

/*----------------------------------------------------------------------------*/

void AliAnalysisTaskGammaPHOSPP::FillTwoPhotonHistograms()
{
  TLorentzVector p1, p2, p12, pv1, pv2, pv12, p11;

  for (Int_t i1 = 0; i1 < fInPHOS-1; i1 ++ ) 
  {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    for (Int_t i2=i1+1; i2<fInPHOS; i2++) 
    {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      p12  = *ph1  + *ph2;
      pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
      Double_t asym  = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));
      Double_t ma12 = p12.M();
      Double_t pt12 = p12.Pt();
      Double_t w = fWeightFunction->Eval(pt12);

      if (ph1->GetNCells()>2 && ph2->GetNCells()>2) 
      {
        FillHistogram("hMassPtA10", ma12 , pt12, w );
        FillHistogram("hMassPtvA10", pv12.M(), pv12.Pt(), w);
        FillHistogram("hMassPtCA10", ma12 , pt12, fEventCentrality + 0.5);
        FillHistogram("hMassSingle_all", ma12,ph1->Pt(), w) ;
        FillHistogram("hMassSingle_all", ma12,ph2->Pt(), w) ;

        FillHistogram("hMinv_all_all",ma12,ph1->Pt()); //!!!!!!!!!!!!!!!!!
	if(!ph1->IsDispOK()) FillHistogram("hMinv_antidisp_all",ma12,ph1->Pt());//!!!!!!!!!!
	if(ph1->IsDispOK()) FillHistogram("hMinv_disp_all",ma12,ph1->Pt());//!!!!!!!!!!

	if(!ph1->IsCPVOK())FillHistogram("hMinv_anticpv_all",ma12,ph1->Pt());//!!!!!!!!!!
	if(ph1->IsCPVOK()) FillHistogram("hMinv_cpv_all",ma12,ph1->Pt());//!!!!!!!!!!

        if(i2 + 1 == fInPHOS )
        {
          FillHistogram("hMinvClustPt_all", ma12, ph1->Pt());
          if(ph1->IsCPVOK())  FillHistogram("hMinvClustPt_cpv",   ma12, ph1->Pt());
          if(ph1->IsDispOK()) FillHistogram("hMinvClustPt_disp",  ma12, ph1->Pt());
          if(ph1->IsCPVOK() && ph1->IsDispOK()) FillHistogram("hMinvClustPt_both",  ma12, ph1->Pt());
        }

      //-----------------------------

	//if(!eventVtxExist)
	  FillHistogram("hMassPtA10nvtx",ma12 ,pt12,w );
	//if(eventVtxExist)
	  FillHistogram("hMassPtA10vtx"  ,ma12 ,pt12,w );
	//if(eventVtxExist && eventV0AND)
	  FillHistogram("hMassPtA10V0AND",ma12 ,pt12,w );
	//if(eventPileup)
	  FillHistogram("hMassPtA10PU"   ,ma12 ,pt12,w );
        if(ph1->IsCPVOK())
          FillHistogram("hMassSingle_cpv",ma12,ph1->Pt(),w) ;
        if(ph2->IsCPVOK())
          FillHistogram("hMassSingle_cpv",ma12,ph2->Pt(),w) ;
        if(ph1->IsDispOK())
          FillHistogram("hMassSingle_disp",ma12,ph1->Pt(),w) ;
        if(ph2->IsDispOK())
          FillHistogram("hMassSingle_disp",ma12,ph2->Pt(),w) ;
        if(ph1->IsCPVOK() && ph1->IsDispOK())
          FillHistogram("hMassSingle_both",ma12,ph1->Pt(),w) ;
        if(ph2->IsCPVOK() && ph2->IsDispOK())
          FillHistogram("hMassSingle_both",ma12,ph2->Pt(),w) ;
        if(ph1->IsCPVOK() && ph2->IsCPVOK())
          FillHistogram("hMassPtCA10_cpv",ma12 ,pt12, fEventCentrality + 0.5);
        if(ph1->IsDispOK() && ph2->IsDispOK())
        {
          FillHistogram("hMassPtCA10_disp",ma12 ,pt12, fEventCentrality + 0.5);
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMassPtCA10_both",ma12 ,pt12, fEventCentrality + 0.5);
        }
        if (asym<0.8) {
          FillHistogram("hMassPtA08",ma12,pt12);
        }
        if (asym<0.7) 
        {
          FillHistogram("hMassPtA07", ma12, pt12);
	  FillHistogram("hMassPtvA07", pv12.M(), pv12.Pt());
          FillHistogram("hMassPtCA07", ma12, pt12, fEventCentrality + 0.5);
	 // if(!eventVtxExist)
	    FillHistogram("hMassPtA07nvtx", ma12, pt12 );
	  //if(eventVtxExist)
	    FillHistogram("hMassPtA07vtx"  , ma12,pt12 );
	  //if(eventVtxExist && eventV0AND)
	    FillHistogram("hMassPtA07V0AND",ma12 ,pt12 );
	 // if(eventPileup)
	    FillHistogram("hMassPtA07PU"   ,ma12 ,pt12 );
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMassPtCA07_cpv",ma12 , pt12, fEventCentrality + 0.5);
          if(ph1->IsDispOK() && ph2->IsDispOK())
          {
            FillHistogram("hMassPtCA07_disp",ma12 ,pt12, fEventCentrality + 0.5);
            if(ph1->IsCPVOK() && ph2->IsCPVOK())
              FillHistogram("hMassPtCA07_both",ma12 ,pt12, fEventCentrality + 0.5);
          }
        }
        if (asym<0.1) 
           FillHistogram("hMassPtA01",ma12,pt12);
	if (TMath::Abs(ma12-0.135)<0.03) 
	   FillHistogram("hAsymPtPi0", asym   ,pt12);
	if (TMath::Abs(ma12-0.547)<0.09) 
	   FillHistogram("hAsymPtEta", asym   ,pt12);

        if (ph1->Module()==1 && ph2->Module()==1) 
        {
	  FillHistogram("hMassPtM1",ma12 ,pt12 );
	  if (TMath::Abs(ma12-0.135)<0.03) FillHistogram("hAsymPtPi0M1",asym   ,pt12);
	}
        if (ph1->Module()==2 && ph2->Module()==2) 
        {
	  FillHistogram("hMassPtM2",ma12 ,pt12 );
	  if (TMath::Abs(ma12-0.135)<0.03) FillHistogram("hAsymPtPi0M2",asym   ,pt12);
	}
        if (ph1->Module()==3 && ph2->Module()==3) 
        {
	  FillHistogram("hMassPtM3",ma12 ,pt12 );
	  if (TMath::Abs(ma12-0.135)<0.03) FillHistogram("hAsymPtPi0M3",asym   ,pt12);
	}
        if ((ph1->Module()==1 && ph2->Module()==2) ||
	    (ph1->Module()==2 && ph2->Module()==1)) 
        {
	  FillHistogram("hMassPtM12",ma12 ,pt12 );
	  if (TMath::Abs(ma12-0.135)<0.03) FillHistogram("hAsymPtPi0M12",asym   ,pt12);
	}
        if ((ph1->Module()==2 && ph2->Module()==3) ||
	    (ph1->Module()==3 && ph2->Module()==2)) 
        {
	  FillHistogram("hMassPtM23",ma12 ,pt12 );
	  if (TMath::Abs(ma12-0.135)<0.03) FillHistogram("hAsymPtPi0M23",asym   ,pt12);
	}
        if ((ph1->Module()==1 && ph2->Module()==3) ||
	    (ph1->Module()==3 && ph2->Module()==1)) FillHistogram("hMassPtM13",ma12 ,pt12 );

	if (TMath::Abs(ph1->EMCz()) < 20. || TMath::Abs(ph2->EMCz()) < 20.)
	  FillHistogram("hMassPt20cm",ma12 ,pt12 );
	if ((TMath::Abs(ph1->EMCz()) > 20. && TMath::Abs(ph1->EMCz()) < 40.) ||
	    (TMath::Abs(ph2->EMCz()) > 20. && TMath::Abs(ph2->EMCz()) < 40.))
	  FillHistogram("hMassPt40cm",ma12 ,pt12 );
	if (TMath::Abs(ph1->EMCz()) > 40. || TMath::Abs(ph2->EMCz()) > 40.)
	  FillHistogram("hMassPt60cm",ma12 ,pt12 );
      }

      if (ph1->GetNCells()>3 && ph2->GetNCells()>3) FillHistogram("hMassPtN3",ma12 ,pt12 );
      if (ph1->GetNCells()>4 && ph2->GetNCells()>4) FillHistogram("hMassPtN4",ma12 ,pt12 );
      if (ph1->GetNCells()>5 && ph2->GetNCells()>5) FillHistogram("hMassPtN5",ma12 ,pt12 );
      if (ph1->GetNCells()>6 && ph2->GetNCells()>6) FillHistogram("hMassPtN6",ma12 ,pt12 );
    } 
  } 
  
}

/*----------------------------------------------------------------------------*/
void AliAnalysisTaskGammaPHOSPP::MixPhotons()
{

  TLorentzVector p1, p2, p12, pv1, pv2, pv12, p11;

  Int_t zvtx = (Int_t)((fVtx5[2] + 10.)/2.) ;
  if(zvtx < 0) zvtx = 0 ;
  if(zvtx > 9) zvtx = 9 ;
  
  Int_t centr = 0;
  
  if(!fPHOSEvents[zvtx][centr]) fPHOSEvents[zvtx][centr]=new TList() ;
  
  TList * prevPHOS = fPHOSEvents[zvtx][centr] ;

  for (Int_t i1=0; i1<fInPHOS; i1++) 
  {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    for(Int_t ev = 0; ev < prevPHOS->GetSize(); ev ++)
    {
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;
      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++)
      {
      AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
      p12  = *ph1  + *ph2;
      pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
      Double_t asym  = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));
      Double_t ma12 = p12.M();
      Double_t pt12 = p12.Pt();
      Double_t wm = fWeightFunction->Eval(pt12);

      if (ph1->GetNCells()>2 && ph2->GetNCells()>2) 
      {
        FillHistogram("hMiMassPtA10",ma12 ,pt12,wm );
        FillHistogram("hMiMassPtvA10",pv12.M(),pv12.Pt(),wm);
        FillHistogram("hMiMassPtCA10",ma12 ,pt12, fEventCentrality+0.5);
        FillHistogram("hMiMassSingle_all",ma12,ph1->Pt(),wm) ;
        FillHistogram("hMiMassSingle_all",ma12,ph2->Pt(),wm) ;
	//if(!eventVtxExist)
	  FillHistogram("hMiMassPtA10nvtx",ma12 ,pt12,wm );
	//if(eventVtxExist)
	  FillHistogram("hMiMassPtA10vtx"  ,ma12 ,pt12,wm );
	//if(eventVtxExist && eventV0AND)
	  FillHistogram("hMiMassPtA10V0AND",ma12 ,pt12,wm );
	//if(eventPileup)
	  FillHistogram("hMiMassPtA10PU"   ,ma12 ,pt12 );

        if(ph1->IsCPVOK())  FillHistogram("hMiMassSingle_cpv", ma12, ph1->Pt(), wm) ;
        if(ph2->IsCPVOK())  FillHistogram("hMiMassSingle_cpv", ma12, ph2->Pt(), wm) ;
        if(ph1->IsDispOK()) FillHistogram("hMiMassSingle_disp", ma12, ph1->Pt(), wm) ;
        if(ph2->IsDispOK()) FillHistogram("hMiMassSingle_disp", ma12, ph2->Pt(), wm) ;
        if(ph1->IsCPVOK() && ph1->IsDispOK()) FillHistogram("hMiMassSingle_both", ma12, ph1->Pt(), wm) ;
        if(ph2->IsCPVOK() && ph2->IsDispOK()) FillHistogram("hMiMassSingle_both", ma12, ph2->Pt(), wm) ;

        if(ph1->IsCPVOK() && ph2->IsCPVOK())
          FillHistogram("hMiMassPtCA10_cpv", ma12 ,pt12, fEventCentrality + 0.5);
        if(ph1->IsDispOK() && ph2->IsDispOK())
        {
          FillHistogram("hMiMassPtCA10_disp", ma12 ,pt12, fEventCentrality + 0.5);
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMiMassPtCA10_both", ma12 ,pt12, fEventCentrality + 0.5);
        }
        if (asym<0.8) 
        {
          FillHistogram("hMiMassPtA08", ma12, pt12);
        }
        if (asym<0.7) 
        {
          FillHistogram("hMiMassPtA07", ma12, pt12);
 	  FillHistogram("hMiMassPtvA07", pv12.M(), pv12.Pt());
	  FillHistogram("hMiMassPtCA07", ma12 , pt12, fEventCentrality + 0.5);
	  //if(!eventVtxExist)
	    FillHistogram("hMiMassPtA07nvtx", ma12, pt12 );
	  //if(eventVtxExist)
	    FillHistogram("hMiMassPtA07vtx"  , ma12, pt12 );
	  //if(eventVtxExist && eventV0AND)
	    FillHistogram("hMiMassPtA07V0AND", ma12, pt12 );
	 // if(eventPileup)
	    FillHistogram("hMiMassPtA07PU"   , ma12, pt12 );
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMiMassPtCA07_cpv", ma12, pt12, fEventCentrality +0.5);
          if(ph1->IsDispOK() && ph2->IsDispOK())
          {
            FillHistogram("hMiMassPtCA07_disp", ma12, pt12, fEventCentrality + 0.5);
            if(ph1->IsCPVOK() && ph2->IsCPVOK())
              FillHistogram("hMiMassPtCA07_both", ma12 , pt12, fEventCentrality + 0.5);
          }
        }
        if (asym<0.1) FillHistogram("hMiMassPtA01",ma12,pt12);
     
        FillHistogram("hMiAsymPt",asym   ,pt12);

        if (ph1->Module()==1 && ph2->Module()==1) FillHistogram("hMiMassPtM1",ma12 ,pt12 );
        if (ph1->Module()==2 && ph2->Module()==2) FillHistogram("hMiMassPtM2",ma12 ,pt12 );
        if (ph1->Module()==3 && ph2->Module()==3) FillHistogram("hMiMassPtM3",ma12 ,pt12 );
        if ((ph1->Module()==1 && ph2->Module()==2) ||
	    (ph1->Module()==2 && ph2->Module()==1)) FillHistogram("hMiMassPtM12",ma12 ,pt12 );
        if ((ph1->Module()==2 && ph2->Module()==3) ||
	    (ph1->Module()==3 && ph2->Module()==2)) FillHistogram("hMiMassPtM23",ma12 ,pt12 );
        if ((ph1->Module()==1 && ph2->Module()==3) ||
	    (ph1->Module()==3 && ph2->Module()==1)) FillHistogram("hMiMassPtM13",ma12 ,pt12 );

	if (TMath::Abs(ph1->EMCz()) < 20. || TMath::Abs(ph2->EMCz()) < 20.)
	  FillHistogram("hMiMassPt20cm",ma12 ,pt12 );
	if ((TMath::Abs(ph1->EMCz()) > 20. && TMath::Abs(ph1->EMCz()) < 40.) ||
	    (TMath::Abs(ph2->EMCz()) > 20. && TMath::Abs(ph2->EMCz()) < 40.))
	  FillHistogram("hMiMassPt40cm",ma12 ,pt12 );
	if (TMath::Abs(ph1->EMCz()) > 40. || TMath::Abs(ph2->EMCz()) > 40.)
	  FillHistogram("hMiMassPt60cm",ma12 ,pt12 );
      }

      if (ph1->GetNCells()>3 && ph2->GetNCells()>3) FillHistogram("hMiMassPtN3",ma12 ,pt12 );
      if (ph1->GetNCells()>4 && ph2->GetNCells()>4) FillHistogram("hMiMassPtN4",ma12 ,pt12 );
      if (ph1->GetNCells()>5 && ph2->GetNCells()>5) FillHistogram("hMiMassPtN5",ma12 ,pt12 );
      if (ph1->GetNCells()>6 && ph2->GetNCells()>6) FillHistogram("hMiMassPtN6",ma12 ,pt12 );
               
      } // end of loop i2
    }
  } 
 
  if(fPHOSEvent->GetEntriesFast()>0)
  {
    prevPHOS->AddFirst(fPHOSEvent) ;
    fPHOSEvent=0;
    if(prevPHOS->GetSize()>100)
    {//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last()) ;
      prevPHOS->RemoveLast() ;
      delete tmp ;
    }
  }
}

//=================================== Returns label ============================

Int_t AliAnalysisTaskGammaPHOSPP::GetPrimaryLabelAtVertex(AliVCluster *clu)
{
   if(!fMCArray) 
     return 0;
      
   Int_t iPrimaryAtVertex = clu->GetLabel();
   AliAODMCParticle *particle0 =  (AliAODMCParticle*) fMCArray->At(iPrimaryAtVertex);
   
   if(particle0-> IsSecondaryFromMaterial())
   {
     Printf("Secondary from the material, Epart = %f, Eclust = %f" , particle0 ->E(), clu->E());
   //  return 0;
   }

   while(TMath::Hypot(particle0 -> Xv(), particle0 -> Yv()) > 1.0)
   {
      FillHistogram("htest",0.5);
      iPrimaryAtVertex = particle0->GetMother();
      particle0 = (AliAODMCParticle*) fMCArray->At(particle0->GetMother());
   }

   Int_t nn = iPrimaryAtVertex;

   if(particle0->GetPdgCode() == 22 || particle0->GetPdgCode() == 11)
   {
     for(Int_t i = 0; i < nn /*fMCArray->GetEntriesFast()*/; i++)
     {
       AliAODMCParticle* particle =  (AliAODMCParticle*) fMCArray->At(i);
       if(particle->GetPdgCode() != 310 && particle->GetPdgCode() != 130) continue;
       Int_t iSecondDaughter = particle->GetDaughterLabel(1); 
       if(iSecondDaughter != iPrimaryAtVertex) continue;
       else
         iPrimaryAtVertex = i;
     }
   }
   else iPrimaryAtVertex = nn;
  
   return iPrimaryAtVertex;   
}

//=============================================================================

Int_t AliAnalysisTaskGammaPHOSPP::GetPrimaryLabel(AliVCluster *clu)
{
   if(!fMCArray) 
     return 0;
      
   return clu->GetLabel();
}

//=================================== TestGamma returns momentum==============

Double_t AliAnalysisTaskGammaPHOSPP::TestGammaPt(AliCaloPhoton *ph)
{

   if(!fMCArray) 
      return 0;

   Int_t iPrimaryAtVertex=ph->GetPrimaryAtVertex();
   AliAODMCParticle *particle0 =  (AliAODMCParticle*) fMCArray->At(iPrimaryAtVertex);

   return (particle0->Pt());

}

//=======================================
Int_t AliAnalysisTaskGammaPHOSPP::TestTrack(AliAODTrack *track)
{

   if(!fMCArray) 
      return 0;

   Int_t TrackLabel=track->GetLabel();

   //if(TrackLabel < -2) return 0;

   AliAODMCParticle *TrackParticle =  (AliAODMCParticle*) fMCArray->At(TrackLabel);
   if(!TrackParticle) 
       return 0;

 //  if (((AliAODMCParticle*) TrackParticle)->IsSecondaryFromWeakDecay()) return 0;
 //  if (!TrackParticle->IsPhysicalPrimary()) return 0;

   Int_t TrackPDG = TrackParticle->GetPdgCode();

   return TrackPDG;
}

//=======================================
Double_t AliAnalysisTaskGammaPHOSPP::Weight(AliAODMCParticle *particleAtVertex)
{
  if(!fMCArray) 
     return 1.0;
   
   if(particleAtVertex-> IsSecondaryFromMaterial())
   {
  //   Printf("Secondary from the material, Epart = %f" , particleAtVertex->E());
     return 0;
   }

   Int_t iPrimaryAtVertex = particleAtVertex->Label();
   
   while(TMath::Hypot(particleAtVertex -> Xv(), particleAtVertex -> Yv()) > 1.0)
   {
      iPrimaryAtVertex = particleAtVertex->GetMother();
      particleAtVertex = (AliAODMCParticle*)fMCArray->At(particleAtVertex->GetMother());
   }

/*
   Int_t nn = iPrimaryAtVertex;

   if(particle->GetPdgCode() == 22 || particle->GetPdgCode() == 11)
   {
     for(Int_t i = 0; i < nn; i++)
     {
       AliAODMCParticle* particle =  (AliAODMCParticle*) fMCArray->At(i);
       if(particle->GetPdgCode() != 310 && particle->GetPdgCode() != 130) continue;
       Int_t iSecondDaughter = particle->GetDaughterLabel(1); 
       if(iSecondDaughter != iPrimaryAtVertex) continue;
       else
         iPrimaryAtVertex = i;
     }
   }
   else 
      iPrimaryAtVertex = nn;
    
        Printf("11111");
   AliAODMCParticle *particleAtVertex = (AliAODMCParticle*)fMCArray->At(iPrimaryAtVertex);      
 */
   if(TMath::Abs(particleAtVertex->GetPdgCode()) == 111)
      fWeightFunction2->SetParameters(0.611073, -0.0222529, 0.190541, -0.416579, 0.396059, 0.611073);
   else  if(TMath::Abs(particleAtVertex->GetPdgCode()) == 221 || 
            TMath::Abs(particleAtVertex->GetPdgCode()) == 331 ||
            TMath::Abs(particleAtVertex->GetPdgCode()) == 223 )
            fWeightFunction2->SetParameters(0.0601459, 0, 4.11665, 0, 6.46838, -0.00319589);  
         else  if(TMath::Abs(particleAtVertex->GetPdgCode()) == 130 ||
                  TMath::Abs(particleAtVertex->GetPdgCode()) == 310 ||
                  TMath::Abs(particleAtVertex->GetPdgCode()) == 311 ||
                  TMath::Abs(particleAtVertex->GetPdgCode()) == 321 )
                  fWeightFunction2->SetParameters(0.708656, 0.355564, -0.00468263, 0.0570132, 0.076876, 0.0382327);
               else if(TMath::Abs(particleAtVertex->GetPdgCode()) > 1000)
                       fWeightFunction2->SetParameters(0.215726, 0.292934, 0.163074, -0.460113, 0.219988, -0.0903996);
                    else fWeightFunction2->SetParameters(1.0, 0., 0., 0., 0., 0.);
                    
   Double_t pt = particleAtVertex->Pt();
   
   return fWeightFunction2->Eval(pt);
} 

//=======================================
Int_t AliAnalysisTaskGammaPHOSPP::TestBC(Double_t tof)
{
  Int_t bc = (Int_t)(TMath::Ceil((tof + fBCgap/2)/fBCgap) - 1);
  return bc;
}

//===============================================
Bool_t AliAnalysisTaskGammaPHOSPP::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  //
  
  //if(!fCheckMCCrossSection) return kTRUE;

  // Fetch the aod also from the input in,
  // have todo it in notify
  
  Float_t xsection = 0;
  Float_t trials   = 1;
  fAvgTrials = -1;
  
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if(!tree) return kFALSE;
  
  TFile *curfile = tree->GetCurrentFile();
  
  if(!curfile) return kFALSE;
  
  if(fCurrFileName == curfile->GetName()) return kFALSE;
  
  fCurrFileName = TString(curfile->GetName());
  
  if(!fh1Xsec||!fh1Trials)
  {
  //  Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
    return kFALSE;
  }
  
  Bool_t ok = PythiaInfoFromFile(fCurrFileName, xsection, trials);
  
  if(!ok) return kFALSE;
  
  fh1Xsec->Fill("<#sigma>",xsection);
  
  // construct a poor man average trials
  Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
  
  if(trials >= nEntries && nEntries > 0.) fAvgTrials = trials/nEntries;
  
  fh1Trials->Fill("#sum{ntrials}",trials);
  
  //printf("AliAnalysisTaskGammaPHOSPP::Notify() - xs %f, trial %f, avg trials %f\n",xsection,trials, fAvgTrials);
  
  if(fDebug) Printf("Reading File %s",fInputHandler->GetTree()->GetCurrentFile()->GetName());
  
  return kTRUE;
}

//_____________________________________________________________________________________________________
Bool_t AliAnalysisTaskGammaPHOSPP::PythiaInfoFromFile(TString file,Float_t & xsec,Float_t & trials)
{
  //
  // get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  // This is to called in Notify and should provide the path to the AOD/ESD file
    
  xsec   = 0;
  trials = 1;
  
  if(file.Contains("root_archive.zip#"))
  {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos  = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  }
  else
  {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  
  //Printf("%s",file.Data());
  
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
  if(!fxsec)
  {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if(!fxsec)
    {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    }
    else
    {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
      if(!key)
      {
        fxsec->Close();
        return kFALSE;
      }
      
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if(!list)
      {
        fxsec->Close();
        return kFALSE;
      }
      
      xsec    = ((TProfile*)list->FindObject("h1Xsec"))  ->GetBinContent(1);
      trials  = ((TH1F*)    list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else
  {
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if(!xtree)
    {
      fxsec->Close();
      return kFALSE;
    }
    
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    trials = ntrials;
    xsec = xsection;
    fxsec->Close();
  }
  
  return kTRUE;
}

//=================

Bool_t AliAnalysisTaskGammaPHOSPP::PhotonWithinPeak(Double_t Minv, Double_t pt)
{
  const Int_t Nbins = 33;

  Double_t xbins[Nbins+1] = {
                             0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2,
                             2.2, 2.4, 2.6, 2.8, 3, 3.2, 
                             3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 7, 8, 10, 12, 14, 16, 
                             18, 20, 25
                            };

  Double_t xmass  [Nbins] = {
                             0.13531, 0.13531,  0.13531,  0.13531,
                             0.13531, 0.135298, 0.135323, 0.135263, 0.135279, 0.135261, 
                             0.135328,0.135201, 0.135349, 0.135275, 0.1353, 0.1354, 
                             0.135014,0.13545,  0.135336, 0.134832, 0.135384, 0.135042, 
                             0.13505, 0.134911, 0.135474, 0.134483, 0.13551, 0.135829, 
                             0.130291,0.134933, 0.145022, 0.139125, 0.135036
                           };

  Double_t xwidth [Nbins] = {
                             5.47384, 5.47384, 5.47384, 5.47384,
                             5.47384, 5.19915, 5.16869, 5.00068, 5.03612, 4.92791, 
		             4.99442, 4.98119, 5.14243, 5.01174, 5.09855, 4.82253, 
                             5.10046, 5.09045, 5.22204, 5.20583, 5.30179, 5.00537, 
                             5.09301, 5.57206, 5.64171, 4.92138, 5.49863, 7.24342, 
                             8.35561, 5.34015, 8.93808, 2.65933, 1.80287
                            };

  Int_t k=0;

  while(pt > xbins[k] && k < Nbins+1)
       k=k+1;
 
  return( Minv < xmass[k] + xwidth[k] && Minv > xmass[k] - xwidth[k]);
}

//=================
void AliAnalysisTaskGammaPHOSPP::TestMatchingTrackPID(AliVCluster *clu1, Double_t pt)
{
  
    Bool_t CPVBit = kFALSE;
    
    CPVBit = clu1->GetEmcCpvDistance() > 2.5;
        
     Bool_t DispBit = clu1->Chi2() < 2.5*2.5;

    Int_t NTracksMatched = clu1->GetNTracksMatched();
    Int_t nmaxMatched;

    Double_t dx = clu1->GetTrackDx(); 
    Double_t dz = clu1->GetTrackDz(); 

    Double_t dist = TMath::Hypot(dx,dz);
    Double_t pBayesMatched[AliPID::kSPECIESC];

    AliPIDCombined *pidcomb=new AliPIDCombined();

    pidcomb->SetDefaultTPCPriors();
    //pidcomb->SetEnablePriors(kFALSE);
    pidcomb->SetSelectedSpecies(AliPID::kSPECIESC);
    pidcomb->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF|AliPIDResponse::kDetITS|AliPIDResponse::kDetTRD);
      
    if(NTracksMatched>0)
    {
        
     Int_t npid;

    FillHistogram("hTracks_matched", 0.5);
    FillHistogram("hDistance",dist);
    AliAODTrack* trackMatched= dynamic_cast<AliAODTrack*>(clu1->GetTrackMatched(0));

    if(trackMatched->TestFilterBit(32) &&  trackMatched->GetTPCsignal() > 0.)
    {

    UInt_t oo=pidcomb->ComputeProbabilities(trackMatched, fPIDResponse, pBayesMatched);

    
    Bool_t pidPion3 = kFALSE , pidKaon3 = kFALSE , pidProton3 = kFALSE , pidElectron3 = kFALSE, pidUndef3 = kFALSE ;
    Bool_t pidPion1 = kFALSE , pidKaon1 = kFALSE , pidProton1 = kFALSE,  pidElectron1 = kFALSE, pidUndef1 = kFALSE;
   
    Float_t  nsigmaElectron =   TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kElectron )) ;
    Float_t  nsigmaPion =       TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kPion )) ;
    Float_t  nsigmaKaon =       TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kKaon ) ) ;
    Float_t  nsigmaProton =     TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kProton ) );

   // smallest sigma

    if(( nsigmaPion < nsigmaKaon) && (nsigmaPion < nsigmaProton) && (nsigmaPion < nsigmaElectron) && (nsigmaPion < 3.)) pidPion3= kTRUE;
    if(( nsigmaProton < nsigmaKaon) && (nsigmaProton < nsigmaPion) && (nsigmaProton < nsigmaElectron) &&  (nsigmaProton < 3.)) pidProton3= kTRUE;
    if(( nsigmaKaon < nsigmaPion) && (nsigmaKaon < nsigmaProton) && (nsigmaKaon < nsigmaElectron) &&  (nsigmaKaon < 3.)) pidKaon3= kTRUE;
    if(( nsigmaElectron < nsigmaKaon) && (nsigmaElectron < nsigmaPion) && (nsigmaElectron < nsigmaProton) && (nsigmaElectron < 3.)) pidElectron3= kTRUE;
    if(!pidPion3 && !pidProton3 && !pidKaon3 && !pidElectron3) pidUndef3=kTRUE;

    if(pidPion3)     FillHistogram("hpid3", 0.5);     
    if(pidProton3)   FillHistogram("hpid3", 1.5);     
    if(pidKaon3)     FillHistogram("hpid3", 2.5);     
    if(pidElectron3) FillHistogram("hpid3", 3.5);     
    if(pidUndef3)    FillHistogram("hpid3", 4.5);     

    if(( nsigmaPion < nsigmaKaon) && (nsigmaPion < nsigmaProton) && (nsigmaPion < nsigmaElectron) && (nsigmaPion < 1.)) pidPion1= kTRUE;
    if(( nsigmaProton < nsigmaKaon) && (nsigmaProton < nsigmaPion) && (nsigmaProton < nsigmaElectron) &&  (nsigmaProton < 1.)) pidProton1= kTRUE;
    if(( nsigmaKaon < nsigmaPion) && (nsigmaKaon < nsigmaProton) && (nsigmaKaon < nsigmaElectron) &&  (nsigmaKaon < 1.)) pidKaon1= kTRUE;
    if(( nsigmaElectron < nsigmaKaon) && (nsigmaElectron < nsigmaPion) && (nsigmaElectron < nsigmaProton) && (nsigmaElectron < 1.)) pidElectron1= kTRUE;
    if(!pidPion1 && !pidProton1 && !pidKaon1 && !pidElectron1) pidUndef1=kTRUE;

    if(pidPion1)     FillHistogram("hpid1", 0.5);     
    if(pidProton1)   FillHistogram("hpid1", 1.5);     
    if(pidKaon1)     FillHistogram("hpid1", 2.5);     
    if(pidElectron1) FillHistogram("hpid1", 3.5);     
    if(pidUndef1)    FillHistogram("hpid1", 4.5);     
    
    nmaxMatched=TMath::LocMax(AliPID::kSPECIESC, pBayesMatched);

    FillHistogram("hTracksofClusts",pt,nmaxMatched+0.5);
    FillHistogram("hEnTrackvsClust",clu1->E(), trackMatched->E());

        
        if(pidPion3 && trackMatched->Charge() > 0)
        {
            FillHistogram("hTracksOfPi_ThreeSigma", pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfPi_ThreeSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit) 
            {
              FillHistogram("hTracksOfPi_ThreeSigma_disp",       pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfPi_ThreeSigma_disp_label", pt, TestTrack(trackMatched) );
            }

        } 
        if(pidPion1 && trackMatched->Charge() > 0)
        {
            FillHistogram("hTracksOfPi_OneSigma", pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfPi_OneSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit) 
            {
              FillHistogram("hTracksOfPi_OneSigma_disp",       pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfPi_OneSigma_disp_label", pt, TestTrack(trackMatched) );
            }
        } 
        if(pidPion3 && trackMatched->Charge() < 0)
        {
            FillHistogram("hTracksOfAntiPi_ThreeSigma", pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfAntiPi_ThreeSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit) 
            {
              FillHistogram("hTracksOfAntiPi_ThreeSigma_disp",       pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfAntiPi_ThreeSigma_disp_label", pt, TestTrack(trackMatched) );
            }
        } 
        if(pidPion1 && trackMatched->Charge() < 0)
        {
            FillHistogram("hTracksOfAntiPi_OneSigma", pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfAntiPi_OneSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit) 
            {
              FillHistogram("hTracksOfAntiPi_OneSigma_disp",       pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfAntiPi_OneSigma_disp_label", pt, TestTrack(trackMatched) );
            }
        } 
        if(pidProton3 && trackMatched->Charge() > 0)
        { 
            FillHistogram("hTracksOfPr_ThreeSigma",pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfPr_ThreeSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit)
            {
              FillHistogram("hTracksOfPr_ThreeSigma_disp_label", pt, TestTrack(trackMatched) );
              FillHistogram("hTracksOfPr_ThreeSigma_disp",pt,clu1->GetEmcCpvDistance());
            }
        } 
        if(pidProton1 && trackMatched->Charge() > 0)
        { 
            FillHistogram("hTracksOfPr_OneSigma",pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfPr_OneSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit)
            {
              FillHistogram("hTracksOfPr_OneSigma_disp_label", pt, TestTrack(trackMatched) );
              FillHistogram("hTracksOfPr_OneSigma_disp",pt,clu1->GetEmcCpvDistance());
            }
        }
        if(pidProton3 && trackMatched->Charge() < 0)
        { 
            FillHistogram("hTracksOfAntiPr_ThreeSigma",pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfAntiPr_ThreeSigma_label",pt,TestTrack(trackMatched));
            if(DispBit) 
            {
              FillHistogram("hTracksOfAntiPr_ThreeSigma_disp",       pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfAntiPr_ThreeSigma_disp_label", pt, TestTrack(trackMatched));
            }
        }
        if(pidProton1 && trackMatched->Charge() < 0)
        { 
            FillHistogram("hTracksOfAntiPr_OneSigma",pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfAntiPr_OneSigma_label",pt,TestTrack(trackMatched));
            if(DispBit) 
            {
              FillHistogram("hTracksOfAntiPr_OneSigma_disp",       pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfAntiPr_OneSigma_disp_label", pt, TestTrack(trackMatched));
            }
        }
        if(pidKaon3 && trackMatched->Charge() > 0)
        { 
            FillHistogram("hTracksOfKa_ThreeSigma",       pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfKa_ThreeSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit) 
            {
              FillHistogram("hTracksOfKa_ThreeSigma_disp",       pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfKa_ThreeSigma_disp_label", pt, TestTrack(trackMatched) );
            }
        } 
        if(pidKaon1 && trackMatched->Charge() > 0)
        { 
            FillHistogram("hTracksOfKa_OneSigma",       pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfKa_OneSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit) 
            {
              FillHistogram("hTracksOfKa_OneSigma_disp",       pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfKa_OneSigma_disp_label", pt, TestTrack(trackMatched) );
            }
        } 
        if(pidKaon3 && trackMatched->Charge() < 0)
        { 
            FillHistogram("hTracksOfAntiKa_ThreeSigma",       pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfAntiKa_ThreeSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit)
            { 
              FillHistogram("hTracksOfAntiKa_ThreeSigma_disp",   pt,clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfAntiKa_ThreeSigma_disp_label", pt, TestTrack(trackMatched) );
            }
        } 
        if(pidKaon1 && trackMatched->Charge() < 0)
        { 
            FillHistogram("hTracksOfAntiKa_OneSigma",       pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfAntiKa_OneSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit)
            { 
              FillHistogram("hTracksOfAntiKa_OneSigma_disp",       pt,clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfAntiKa_OneSigma_disp_label", pt, TestTrack(trackMatched) );
            }
        } 
        if(pidElectron3 && trackMatched->Charge() > 0)
        {
            FillHistogram("hTracksOfAntiBeta_ThreeSigma",	  pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfAntiBeta_ThreeSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit)
            {
              FillHistogram("hTracksOfAntiBeta_ThreeSigma_disp",	 pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfAntiBeta_ThreeSigma_disp_label", pt, TestTrack(trackMatched) );
            }
	}
	if(pidElectron1 && trackMatched->Charge() > 0)
        {
            FillHistogram("hTracksOfAntiBeta_OneSigma",	pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfAntiBeta_OneSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit)
            {
              FillHistogram("hTracksOfAntiBeta_OneSigma_disp",       pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfAntiBeta_OneSigma_disp_label", pt, TestTrack(trackMatched) );
            }
	}
        if(pidElectron3 && trackMatched->Charge() < 0)
        {
            FillHistogram("hTracksOfBeta_ThreeSigma",	  pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfBeta_ThreeSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit)
            {
              FillHistogram("hTracksOfBeta_ThreeSigma_disp",	 pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfBeta_ThreeSigma_disp_label", pt, TestTrack(trackMatched) );
            }
	}
	if(pidElectron1 && trackMatched->Charge() < 0)
        {
            FillHistogram("hTracksOfBeta_OneSigma",	pt,clu1->GetEmcCpvDistance());
            FillHistogram("hTracksOfBeta_OneSigma_label", pt, TestTrack(trackMatched) );
            if(DispBit)
            {
              FillHistogram("hTracksOfBeta_OneSigma_disp",       pt, clu1->GetEmcCpvDistance());
              FillHistogram("hTracksOfBeta_OneSigma_disp_label", pt, TestTrack(trackMatched) );
            }
	}
        {
            FillHistogram("hTracksOfOthers",pt,clu1->GetEmcCpvDistance());
            if(DispBit) 
              FillHistogram("hTracksOfOthers_disp",pt,clu1->GetEmcCpvDistance());
        }

        
        if(CPVBit)
        {
            if(nmaxMatched==2)
            { 
                FillHistogram("hTracksOfPiClose",pt);
                if(DispBit) 
                    FillHistogram("hTracksOfPiCloseDispOK",pt);
            }
            else
            if(nmaxMatched==4)
            { 
                FillHistogram("hTracksOfPrClose",pt);
                if(DispBit) 
                    FillHistogram("hTracksOfPrCloseDispOK",pt);
            }
            else
            if(nmaxMatched==3)
            { 
                FillHistogram("hTracksOfKaClose",pt);
                if(DispBit) 
                    FillHistogram("hTracksOfKaCloseDispOK",pt);
            }
            else
            {
                FillHistogram("hTracksOfOthersClose",pt);
                FillHistogram("hTracksOfOthersCloseDispOK",pt);
            }               
         }
      }
   }
}

/*----------------------------------------------------------------------------*/

void AliAnalysisTaskGammaPHOSPP::PHOSvsEMCALClusters()
{
   Int_t multPHOS = 0 , multEMCAL = 0;

   AliAODCaloCluster *clu2; 
   for(Int_t ic = 0; ic < fEvent->GetNumberOfCaloClusters(); ic++)
   {
      clu2 = fEvent->GetCaloCluster(ic);
      if(clu2->IsPHOS() ) multPHOS  =  multPHOS  + 1;
      if(clu2->IsEMCAL()) multEMCAL =  multEMCAL + 1;
   }

   Printf("There are %d caloclusters in this event: %d in PHOS, %d in EMCAL", 
   fEvent->GetNumberOfCaloClusters(), multPHOS, multEMCAL);

   FillHistogram("hPHOSvsEMCAL", 0.5, multPHOS );
   FillHistogram("hPHOSvsEMCAL", 1.5, multEMCAL);
}
