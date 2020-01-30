#if !defined (__CINT__) || defined (__CLING__)
#include <string>
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskGammaPHOSPP.h"
#include "AliCaloPhoton.h"
#include "AliTender.h"
#endif
void LoadEnv();
void runTaskGammaPHOS13TeV( Bool_t isMC    =  kFALSE,
                           TString period  = "LHC16g", 
                           TString runmode = "terminate",
                           Bool_t  local   = kTRUE,
                           TString pass    = "pass1" )
{
      
    LoadEnv();
    Bool_t gridTest = kFALSE;
    TGrid::Connect("alien://");

    AliAnalysisManager *mgr = new AliAnalysisManager("GammaAnalysis");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
    
    const char*  tenderOption = isMC ? "Run2NoNcellCut" : "Run2Tune";
    Int_t        tenderPass   = 1;
    Int_t        recoPass   = 1;

    AliPHOSTenderTask *tenderPHOS = reinterpret_cast<AliPHOSTenderTask *>(gInterpreter->ExecuteMacro(Form("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C(\"%s\", \"%s\", \"%s\", %d, %d)", "PHOSTenderTask","PHOStender", tenderOption, tenderPass, isMC)));

    TString nonlinearity = isMC ? "Run2TuneMCNoNcell" : "Run2Tune";
    tenderPHOS->GetPHOSTenderSupply()->SetNonlinearityVersion(nonlinearity) ;
 
    if(isMC)
      tenderPHOS->GetPHOSTenderSupply()->ApplyZeroSuppression(0.020);


    TMacro addresp(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
    addresp.Exec(Form("%d", isMC));
    addresp.Exec(Form("%d", recoPass));

    gROOT->LoadMacro("AliCaloPhoton.cxx++g");
    gROOT->LoadMacro("AliAnalysisTaskGammaPHOSPP.cxx++g");
    AliAnalysisTaskGammaPHOSPP *task = 
       reinterpret_cast<AliAnalysisTaskGammaPHOSPP*>((gInterpreter->ExecuteMacro(
        Form("AddTaskGammaPHOSPP.C(%d)", isMC))));

     TMacro physseladd(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"));
     AliPhysicsSelectionTask *physseltask = reinterpret_cast<AliPhysicsSelectionTask *>(physseladd.Exec(Form("%d", isMC)));
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    // add your task to the manager


    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) 
    {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        if(isMC)
        {
           chain->Add("alien:///alice/sim/2017/LHC17d20a1_extra/258270/AOD/045/AliAOD.root");        
           chain->Add("alien:///alice/sim/2017/LHC17d20a1_extra/258270/AOD/046/AliAOD.root");
        }   
        else
           chain->Add("alien:///alice/data/2016/LHC16g/000254128/pass1/AOD208/0031/AliAOD.root");

        mgr->StartAnalysis("local", chain);
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskGammaPHOSPP.h AliAnalysisTaskGammaPHOSPP.cxx AliCaloPhoton.h AliCaloPhoton.cxx  libTender.so libTenderSupplies.so libPWGGAPHOSTasks.so");
        alienHandler->SetAnalysisSource("AliAnalysisTaskGammaPHOSPP.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20190821_ROOT6-1");
        alienHandler->SetAPIVersion("V1.1x");
        
        
        // select the input data        
        Int_t year = 2015;
        if(period.Contains("LHC15")) year = 2015;
        if(period.Contains("LHC16")) year = 2016;
        if(period.Contains("LHC17")) year = 2017;
        if(period.Contains("LHC18")) year = 2018; 
                             
        if(isMC)
        {
           alienHandler->SetGridDataDir(Form("/alice/sim/%d/%s", year, period.Data()));
           alienHandler->SetDataPattern("AOD/*/AliAOD.root");
        }
        else
        {
           alienHandler->SetGridDataDir(Form("/alice/data/%d/%s", year, period.Data()));
          // alienHandler->SetDataPattern(Form("/%s/AOD208/*/AliAOD.root", pass.Data()));
           alienHandler->SetDataPattern(Form("/%s/AOD/*/AliAOD.root", pass.Data()));
           alienHandler->SetRunPrefix("000");
        }
        // define the output folders
        if(isMC)
          alienHandler->SetGridWorkingDir(Form("pp_Analysis/%s", period.Data()));
        else
          alienHandler->SetGridWorkingDir(Form("pp_Analysis/%s_%s", period.Data(), pass.Data()));            
        alienHandler->SetGridOutputDir("output");

        // Add runs by numbers
        Int_t Nrun[500], nn = 0;
        Int_t evN;
        ifstream ff;
        if(isMC)
          ff.open(Form("datasets/runs_%s.dat", period.Data()));
        else  
          ff.open(Form("datasets/%s-%s.txt", period.Data(), pass.Data()));
        while( !ff.eof() )
        {
         ff>>Nrun[nn];  
         nn = nn + 1;
        }
        ff.close();
 
        for(Int_t i = 0; i < nn; i++)       
          alienHandler->AddRunNumber(Nrun[i]);
       
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("myTask.sh");
        alienHandler->SetJDLName("myTask.jdl");        
        alienHandler->SetTTL(18000);

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(kTRUE);

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) 
        {
            alienHandler->SetNtestFiles(1);
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
         } 
         else 
         {
            alienHandler->SetRunMode(runmode.Data());
            mgr->StartAnalysis("grid");
         }
    }
}
void LoadEnv()
{

  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  
  //load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPHOSUtils.so"); 
  gSystem->Load("libPWGGAUtils.so");

   //PHOS Tender

  gSystem->Load("libTender.so");
  gSystem->Load("libTenderSupplies.so");
  gSystem->Load("libPWGGAPHOSTasks.so"); 
  
  #if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
    gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/PHOS");
  #else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/PHOS");
  #endif


}
