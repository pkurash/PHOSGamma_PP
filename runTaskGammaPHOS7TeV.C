#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskGammaPHOSPP.h"
#include "AliCaloPhoton.h"
#include "AliTender.h"
#include "TApplication.h"
#endif
void LoadEnv(); 
void runTaskGammaPHOS7TeV( Bool_t isMC    =  kFALSE,
                           TString period  = "LHC10b", 
                           TString runmode = "terminate",
                           Bool_t  local    = kTRUE)
{
    LoadEnv();
    TGrid::Connect("alien://");  
    Bool_t gridTest = kFALSE;

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("GammaAnalysis");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
    
    const char*  tenderOption = isMC ? "LHC14j4b" : "";
    Int_t        tenderPass   = isMC ? 1 : 4;
    Int_t        recoPass     = 4;

     AliPHOSTenderTask *tender = reinterpret_cast<AliPHOSTenderTask *>(
        gInterpreter->ExecuteMacro(
        Form("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C(\"%s\",\"%s\",\"%s\",%d, %d)", 
             "PHOSTenderTask","PHOStender", tenderOption, tenderPass, isMC)));

    AliPHOSTenderSupply *supply = tender->GetPHOSTenderSupply();
    supply->ForceUsingBadMap("alien:///alice/cern.ch/user/p/pkurash/BadMap_LHC10ef_Majority300416.root");
    
    TString nonlinearity = isMC ? "MC" : "Default";
    supply->SetNonlinearityVersion(nonlinearity); 
    if(isMC)
    {
       Double_t NonlinPar[3]={1.008, 0.045, 0.4};
       supply->SetNonlinearityParams(3, NonlinPar);
    }

    // Use custom Zero Suppression threshold if needed     
    if(isMC)
    {
      Double_t zs_threshold = 0.020;
      supply->ApplyZeroSuppression(zs_threshold); 
    }

    TMacro addresp(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
    addresp.Exec(Form("%d", isMC));
    addresp.Exec(Form("%d", recoPass));

    gROOT->LoadMacro("AliCaloPhoton.cxx++g");
    gROOT->LoadMacro(" AliAnalysisTaskGammaPHOSPP.cxx++g");
    AliAnalysisTaskGammaPHOSPP *task = 
       reinterpret_cast<AliAnalysisTaskGammaPHOSPP*>((gInterpreter->ExecuteMacro(
        Form("AddTaskGammaPHOSPP.C(%d)", isMC))));

     TMacro physseladd(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"));
     AliPhysicsSelectionTask *physseltask = reinterpret_cast<AliPhysicsSelectionTask *>(physseladd.Exec(Form("%d", isMC)));
    task->SelectCollisionCandidates(AliVEvent::kMB);      
  
        // Apply recalibration
    if(!isMC)
    {
       task->SetRecalib(1, 0.9822696);
       task->SetRecalib(2, 0.9861288);
       task->SetRecalib(3, 1.0072);
    }     


     if(!mgr->InitAnalysis()) return;
     mgr->SetDebugLevel(2);
     mgr->PrintStatus();
     mgr->SetUseProgressBar(1, 25);


    if(local) 
    {
        TChain* chain = new TChain("aodTree");
        if(isMC)
          chain->Add("alien:///alice/sim/2014/LHC14j4e/129647/AOD/074/AliAOD.root");
        else
          chain->Add("alien:///alice/data/2010/LHC10d/000126285/pass4/AOD172/0001/AliAOD.root");

        mgr->StartAnalysis("local", chain);
    } else 
      {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        alienHandler->SetAdditionalLibs("AliAnalysisTaskGammaPHOSPP.h AliAnalysisTaskGammaPHOSPP.cxx AliCaloPhoton.h AliCaloPhoton.cxx  libTender.so libTenderSupplies.so libPWGGAPHOSTasks.so");
        alienHandler->SetAnalysisSource("AliAnalysisTaskGammaPHOSPP.cxx");
        
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20190322_ROOT6-1");
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        
        // select the input data        
        Int_t year = 2010;
        if(period.Contains("LHC10")) year = 2010;
        if(period.Contains("LHC11")) year = 2011;
        if(period.Contains("LHC12")) year = 2012;
        if(period.Contains("LHC13")) year = 2013; 
        if(period.Contains("LHC14")) year = 2014;
        if(period.Contains("LHC15")) year = 2015;   
        cout << "year = " << year << ", period = " << period.Data() << endl;                                                         
        if(isMC)
        {
           alienHandler->SetGridDataDir(Form("/alice/sim/%d/%s", year, period.Data()));
           alienHandler->SetDataPattern("AOD/*/AliAOD.root");
        }
        else
        {
           alienHandler->SetGridDataDir(Form("/alice/data/%d/%s", year, period.Data()));
           alienHandler->SetDataPattern("pass4/AOD172/*/AliAOD.root");
           alienHandler->SetRunPrefix("000");
        }        
                
        // define the output folders
        alienHandler->SetGridWorkingDir(Form("pp_Analysis/%s", period.Data()));
        alienHandler->SetGridOutputDir("output");
        
        //define run numbers
        Int_t evN[500], nn = 0;
        ifstream ff;
        ff.open(Form("datasets/runs_%s.list", period.Data()));
        while( !ff.eof() )
        {
         ff>>evN[nn];
         nn = nn + 1;
        }
        ff.close();
        for(Int_t  i = 0; i < nn; i ++)
            alienHandler->AddRunNumber(evN[i]);

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("myTask.sh");
        alienHandler->SetJDLName("myTask.jdl");
         alienHandler->SetTTL(21600);

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(kTRUE);

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) 
        {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } 
           else 
           {
              // else launch the full grid analysis
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
  
      // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
    gInterpreter->ProcessLine(".include $ALICE_PHYSCIS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/PHOS");

#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/PHOS");
#endif
     
}
