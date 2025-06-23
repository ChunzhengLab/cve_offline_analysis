// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "TInterpreter.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "TMacro.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliMultSelectionTask.h"
#include "TSystem.h"
#include "TChain.h"
#include "AliAnalysisAlien.h"
#include "AliPhysicsSelectionTask.h"
#include "AliAnalysisTaskMyTaskMyMCLambda.h"
#include <algorithm>
#include <iostream>

void runAnalysis()
{
  // Task Name
  TString taskName = "EfficiencyCorrection2025";
  // Dataset
  TString dataset = "LHC18q";
  bool isMC = 1;

  // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
  // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
  bool local      = 0;
  bool gridTest   = 0;
  bool gridRun    = 0;
  bool gridMerge  = 0;
  bool localMerge = 1;

  std::vector<bool> boolVector = {local, gridTest, gridRun, gridMerge, localMerge};
  int trueCount = std::count(boolVector.begin(), boolVector.end(), true);
  if (trueCount != 1)
  {
    std::cout << "开关只能开一个!" << std::endl;
    return;
  }
  if (gridMerge || localMerge)
    gridRun = kFALSE;
  if (gridRun)
  {
    gridTest = kFALSE;
    local = kFALSE;
  }

  Int_t groupNumber = 1;

  // since we will compile a class, tell root where to look for headers
#if !defined(__CINT__) || defined(__CLING__)
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif

  // create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
  //AliESDInputHandler *esdH = new AliESDInputHandler();
  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  //mgr->SetInputEventHandler(esdH);

  TMacro PIDadd(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
  AliAnalysisTaskPIDResponse *PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse *>(PIDadd.Exec());

  if (dataset.EqualTo("LHC15o") || dataset.EqualTo("LHC18q") || dataset.EqualTo("LHC18r"))
  {
    // TMacro physicsSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"));
    // AliPhysicsSelectionTask* physicsSelectionTask = reinterpret_cast<AliPhysicsSelectionTask*>(physicsSelection.Exec());
    TMacro multSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
    AliMultSelectionTask *multSelectionTask = reinterpret_cast<AliMultSelectionTask *>(multSelection.Exec());
  }

  // compile the class and load the add task macro
  // here we have to differentiate between using the just-in-time compiler
  // from root6, or the interpreter of root5
#if !defined(__CINT__) || defined(__CLING__)
  gInterpreter->LoadMacro("AliAnalysisTaskMyTaskMyMCLambda.cxx++");
  AliAnalysisTaskMyTaskMyMCLambda *task = nullptr;
  task = reinterpret_cast<AliAnalysisTaskMyTaskMyMCLambda*>(gInterpreter->ExecuteMacro("AddTaskMyTaskMyMCLambda.C"));
#else
  gROOT->LoadMacro("AliAnalysisTaskMyTaskMyMCLambda.cxx++");
  gROOT->LoadMacro("AddTaskMyTaskMyMCLambda.C");
  AliAnalysisTaskMyTaskMyMCLambda *task = AddTaskMyTaskMyMCLambda();
#endif

  if (!mgr->InitAnalysis())
    return;
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  mgr->SetUseProgressBar(1, 25);

  if (local)
  {
    // if you want to run locally, we need to define some input
    TChain *chain = new TChain("aodTree");
    // add a few files to the chain (change this so that your local files are added)
    // if (dataset.EqualTo("LHC10h"))
    //   chain->Add("/Users/wangchunzheng/alice/data/2010/LHC10h/000139510/ESDs/pass2/AOD160/0247/AliAOD.root");
    // if (dataset.EqualTo("LHC15o"))
    //   chain->Add("/Users/wangchunzheng/alice/data/2015/LHC15o/000245151/pass2/AOD252/0008/AliAOD.root");
    // if (dataset.EqualTo("LHC18q"))
    //   chain->Add("/Users/wangchunzheng/alice/data/2018/LHC18q/000295588/pass3/AOD252/AOD/001/AliAOD.root");
    // if (dataset.EqualTo("LHC18r"))
    //   chain->Add("/Users/wangchunzheng/alice/data/2018/LHC18r/000296691/pass3/AOD252/0001/AliAOD.root");
    if(isMC) {
        if (dataset.EqualTo("LHC18q") || dataset.EqualTo("LHC18r"))
            chain->Add("/Users/wangchunzheng/alice/sim/2020/LHC20e3a/295585/AOD243/0001/AliAOD.root");
    }
    else {
        if (dataset.EqualTo("LHC18q"))
            chain->Add("/Users/wangchunzheng/alice/data/2018/LHC18q/000295588/pass3/AOD252/AOD/001/AliAOD.root");
        if (dataset.EqualTo("LHC18r"))
            chain->Add("/Users/wangchunzheng/alice/data/2018/LHC18r/000296691/pass3/AOD252/0001/AliAOD.root");
    }
    // start the analysis locally, reading the events from the tchain
    mgr->StartAnalysis("local", chain);
  }
  else
  {
    // if we want to run on grid, we create and configure the plugin
    AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
    // also specify the include (header) paths on grid
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    // make sure your source files get copied to grid
    alienHandler->SetAdditionalLibs("AliAnalysisTaskMyTaskMyMCLambda.cxx AliAnalysisTaskMyTaskMyMCLambda.h");
    alienHandler->SetAnalysisSource("AliAnalysisTaskMyTaskMyMCLambda.cxx");
    // select the aliphysics version. all other packages
    // are LOADED AUTOMATICALLY!
    alienHandler->SetAliPhysicsVersion("vAN-20220316_ROOT6-1");
    // set the Alien API version
    alienHandler->SetAPIVersion("V1.1x");

    if (dataset.EqualTo("LHC10h"))
    {
      // 10h
      //  select the input data
      alienHandler->SetGridDataDir("/alice/data/2010/LHC10h");
      alienHandler->SetDataPattern("ESDs/pass2/AOD160/*/AliAOD.root");
      // MC has no prefix, data has prefix 000
      alienHandler->SetRunPrefix("000");
      // runnumber
      // 1
      if (groupNumber == 1)
      {
        alienHandler->AddRunNumber(139510);
        alienHandler->AddRunNumber(139507);
        alienHandler->AddRunNumber(139505);
        alienHandler->AddRunNumber(139503);
        alienHandler->AddRunNumber(139465);
      }
    }
    if (dataset.EqualTo("LHC15o"))
    {
      // select the input data
      if (isMC)
      {
        // /alice/sim/2016/LHC16g1/246994/AOD198/0162
        // MC has no prefix, data has prefix 000
        alienHandler->SetGridDataDir("/alice/sim/2016/LHC16g1");
        alienHandler->SetDataPattern("AOD198/*/AliAOD.root");
        alienHandler->SetRunPrefix("");
      } else {
        alienHandler->SetGridDataDir("/alice/data/2015/LHC15o");
        alienHandler->SetDataPattern("pass2/AOD252/AOD/*/AliAOD.root");
        alienHandler->SetRunPrefix("000");
      }
      if (groupNumber == 1)
      {
        // "246994", "246991", "246989", "246984", "246982", "246948", "246945", "246928", "246871", "246870",
        // "246867", "246865", "246864", "246859", "246858", "246851", "246847", "246846", "246845", "246844",
        // "246810", "246809", "246808", "246807", "246805", "246804", "246766", "246765", "246763", "246760",
        // "246759", "246758", "246757", "246751", "246750", "246434", "246431", "246424", "246392", "246391",
        // "246276", "246275", "246272", "246271", "246225", "246222", "246217", "246185", "246182", "246181",
        // "246180", "246178", "246153", "246152", "246151", "246148", "246115", "246113", "246089", "246087",
        // "246053", "246052", "246049", "246048", "246042", "246037", "246036", "246012", "246003", "246001",
        // "245963", "245954", "245952", "245949", "245923", "245833", "245831", "245829", "245793", "245785",
        // "245775", "245766", "245759", "245752", "245731", "245729", "245705", "245702", "245692", "245683",
        // "245554", "245545", "245544", "245543", "245542", "245540", "245535", "245507", "245505", "245504",
        // "245501", "245497", "245496", "245454", "245453", "245450", "245446", "245441", "245411", "245410",
        // "245409", "245407", "245401", "245397", "245396", "245353", "245349", "245347", "245346", "245345",
        // "245343", "245259", "245233", "245232", "245231", "245152", "245151", "245146", "245145", "245068",
        // "245066", "245064", "244983", "244982", "244980", "244975", "244918", "244917"

        //随机选20个run number
        alienHandler->AddRunNumber(246276);
        alienHandler->AddRunNumber(246180);
        alienHandler->AddRunNumber(246053);
        alienHandler->AddRunNumber(245963);
        alienHandler->AddRunNumber(245775);
        alienHandler->AddRunNumber(245545);
        alienHandler->AddRunNumber(245233);
        alienHandler->AddRunNumber(245145);
        alienHandler->AddRunNumber(244982);
        alienHandler->AddRunNumber(244917);

        alienHandler->AddRunNumber(246994);
        alienHandler->AddRunNumber(246991);
        alienHandler->AddRunNumber(246989);
        alienHandler->AddRunNumber(246984);
        alienHandler->AddRunNumber(246982);
        alienHandler->AddRunNumber(246948);
        alienHandler->AddRunNumber(246945);
        alienHandler->AddRunNumber(246928);
        alienHandler->AddRunNumber(246871);
        alienHandler->AddRunNumber(246870);
      }
    }
    if (dataset.EqualTo("LHC18q"))
    {
      // select the input data
      if(isMC) {
        alienHandler->SetGridDataDir("/alice/sim/2020/LHC20e3a");
        alienHandler->SetDataPattern("AOD243/*/AliAOD.root");
        alienHandler->SetRunPrefix("");
      } else {
        alienHandler->SetGridDataDir("/alice/data/2018/LHC18q");
        alienHandler->SetDataPattern("pass3/AOD252/AOD/*/AliAOD.root");
        alienHandler->SetRunPrefix("000");
      }
      // MC has no prefix, data has prefix 000
      if (groupNumber == 1)
      {
         // "296623","296622","296621","296619","296618","296616","296615","296594","296553","296552",
         // "296551","296550","296548","296547","296516","296512","296511","296510","296509","296472",
         // "296433","296424","296423","296420","296419","296415","296414","296383","296381","296380",
         // "296379","296378","296377","296376","296375","296312","296309","296304","296303","296280",
         // "296279","296273","296270","296269","296247","296246","296244","296243","296242","296241",
         // "296240","296198","296197","296196","296195","296194","296192","296191","296143","296142",
         // "296135","296134","296133","296132","296123","296074","296066","296065","296063","296062",
         // "296060","296016","295942","295941","295937","295936","295913","295910","295909","295861",
         // "295860","295859","295856","295855","295854","295853","295831","295829","295826","295825",
         // "295822","295819","295818","295816","295791","295788","295786","295763","295762","295759",
         // "295758","295755","295754","295725","295723","295721","295719","295718","295717","295714",
         // "295712","295676","295675","295673","295668","295667","295666","295615","295612","295611",
         // "295610","295589","295588","295586","295585"

        //随机选50个run number,一次性选完
        alienHandler->AddRunNumber(296623);
        alienHandler->AddRunNumber(296622);
        alienHandler->AddRunNumber(296621);
        alienHandler->AddRunNumber(296619);
        alienHandler->AddRunNumber(296618);
        alienHandler->AddRunNumber(296616);
        alienHandler->AddRunNumber(296615);
        alienHandler->AddRunNumber(296594);
        alienHandler->AddRunNumber(296553);
        alienHandler->AddRunNumber(296552);
        alienHandler->AddRunNumber(296551);
        alienHandler->AddRunNumber(296550);
        alienHandler->AddRunNumber(296548);
        alienHandler->AddRunNumber(296547);
        alienHandler->AddRunNumber(296516);
        alienHandler->AddRunNumber(296512);
        alienHandler->AddRunNumber(296511);
        alienHandler->AddRunNumber(296510);
        alienHandler->AddRunNumber(296509);
        alienHandler->AddRunNumber(296472);
        alienHandler->AddRunNumber(296433);
        alienHandler->AddRunNumber(296424);
        alienHandler->AddRunNumber(296423);
        alienHandler->AddRunNumber(296420);
        alienHandler->AddRunNumber(296419);
        alienHandler->AddRunNumber(296415);
        alienHandler->AddRunNumber(296414);
        alienHandler->AddRunNumber(296383);
        alienHandler->AddRunNumber(296381);
        alienHandler->AddRunNumber(296380);
        alienHandler->AddRunNumber(296379);
        alienHandler->AddRunNumber(296378);
        alienHandler->AddRunNumber(296377);
        alienHandler->AddRunNumber(296376);
        alienHandler->AddRunNumber(296375);
        alienHandler->AddRunNumber(296312);
        alienHandler->AddRunNumber(296309);
        alienHandler->AddRunNumber(296304);
        alienHandler->AddRunNumber(296303);
        alienHandler->AddRunNumber(296280);
        alienHandler->AddRunNumber(296279);
        alienHandler->AddRunNumber(296273);
        alienHandler->AddRunNumber(296270);
        alienHandler->AddRunNumber(296269);
        alienHandler->AddRunNumber(296247);
        alienHandler->AddRunNumber(296246);
        alienHandler->AddRunNumber(296244);
        alienHandler->AddRunNumber(296243);
        alienHandler->AddRunNumber(296242);
        alienHandler->AddRunNumber(296241);
        alienHandler->AddRunNumber(296240);
        alienHandler->AddRunNumber(296198);
        alienHandler->AddRunNumber(296197);
        alienHandler->AddRunNumber(296196);
        alienHandler->AddRunNumber(296195);
        alienHandler->AddRunNumber(296194);
        alienHandler->AddRunNumber(296192);
        alienHandler->AddRunNumber(296191);
        alienHandler->AddRunNumber(296143);
        alienHandler->AddRunNumber(296142);

        // alienHandler->AddRunNumber(296135);
        // alienHandler->AddRunNumber(296134);
        // alienHandler->AddRunNumber(296133);
        // alienHandler->AddRunNumber(296132);
        // alienHandler->AddRunNumber(296123);
        // alienHandler->AddRunNumber(296074);
        // alienHandler->AddRunNumber(296066);
        // alienHandler->AddRunNumber(296065);
        // alienHandler->AddRunNumber(296063);
        // alienHandler->AddRunNumber(296062);
        // alienHandler->AddRunNumber(296060);
        // alienHandler->AddRunNumber(296016);
        // alienHandler->AddRunNumber(295942);
        // alienHandler->AddRunNumber(295941);
        // alienHandler->AddRunNumber(295937);
        // alienHandler->AddRunNumber(295936);
        // alienHandler->AddRunNumber(295913);
        // alienHandler->AddRunNumber(295910);
        // alienHandler->AddRunNumber(295909);
        // alienHandler->AddRunNumber(295861);
        // alienHandler->AddRunNumber(295860);
        // alienHandler->AddRunNumber(295859);
        // alienHandler->AddRunNumber(295856);
        // alienHandler->AddRunNumber(295855);
        // alienHandler->AddRunNumber(295854);
        // alienHandler->AddRunNumber(295853);
        // alienHandler->AddRunNumber(295831);
        // alienHandler->AddRunNumber(295829);
        // alienHandler->AddRunNumber(295826);
        // alienHandler->AddRunNumber(295825);
        // alienHandler->AddRunNumber(295822);
        // alienHandler->AddRunNumber(295819);
        // alienHandler->AddRunNumber(295818);
        // alienHandler->AddRunNumber(295816);
        // alienHandler->AddRunNumber(295791);
        // alienHandler->AddRunNumber(295788);
        // alienHandler->AddRunNumber(295786);
        // alienHandler->AddRunNumber(295763);
        // alienHandler->AddRunNumber(295762);
        // alienHandler->AddRunNumber(295759);
        // alienHandler->AddRunNumber(295758);
        // alienHandler->AddRunNumber(295755);
        // alienHandler->AddRunNumber(295754);
        // alienHandler->AddRunNumber(295725);
        // alienHandler->AddRunNumber(295723);
        // alienHandler->AddRunNumber(295721);
        // alienHandler->AddRunNumber(295719);
        // alienHandler->AddRunNumber(295718);
        // alienHandler->AddRunNumber(295717);
        // alienHandler->AddRunNumber(295714);
        // alienHandler->AddRunNumber(295712);
        // alienHandler->AddRunNumber(295676);
        // alienHandler->AddRunNumber(295675);
        // alienHandler->AddRunNumber(295673);
        // alienHandler->AddRunNumber(295668);
        // alienHandler->AddRunNumber(295667);
        // alienHandler->AddRunNumber(295666);
        // alienHandler->AddRunNumber(295615);
        // alienHandler->AddRunNumber(295612);
        // alienHandler->AddRunNumber(295611);
        // alienHandler->AddRunNumber(295610);
        // alienHandler->AddRunNumber(295589);
        // alienHandler->AddRunNumber(295588);
        // alienHandler->AddRunNumber(295586);
        // alienHandler->AddRunNumber(295585);
      }
    }
    if (dataset.EqualTo("LHC18r"))
    {
      // select the input data
      if(isMC) {
        alienHandler->SetGridDataDir("/alice/sim/2020/LHC20e3a");
        alienHandler->SetDataPattern("AOD243/*/AliAOD.root");
        alienHandler->SetRunPrefix("");
      } else {
        alienHandler->SetGridDataDir("/alice/data/2018/LHC18r");
        alienHandler->SetDataPattern("pass3/AOD252/AOD/*/AliAOD.root");
        alienHandler->SetRunPrefix("000");
      }
      if (groupNumber == 1)
      {
        alienHandler->AddRunNumber(297595);
        alienHandler->AddRunNumber(297590);
        alienHandler->AddRunNumber(297588);
        alienHandler->AddRunNumber(297558);
        alienHandler->AddRunNumber(297544);
        alienHandler->AddRunNumber(297542);
        alienHandler->AddRunNumber(297541);
        alienHandler->AddRunNumber(297540);
        alienHandler->AddRunNumber(297537);
        alienHandler->AddRunNumber(297512);
        alienHandler->AddRunNumber(297483);
        alienHandler->AddRunNumber(297479);
        alienHandler->AddRunNumber(297452);
        alienHandler->AddRunNumber(297451);
        alienHandler->AddRunNumber(297450);
        alienHandler->AddRunNumber(297446);
        alienHandler->AddRunNumber(297442);
        alienHandler->AddRunNumber(297441);
        alienHandler->AddRunNumber(297415);
        alienHandler->AddRunNumber(297414);
        alienHandler->AddRunNumber(297413);
        alienHandler->AddRunNumber(297406);
        alienHandler->AddRunNumber(297405);
        alienHandler->AddRunNumber(297380);
        alienHandler->AddRunNumber(297379);
        alienHandler->AddRunNumber(297372);
        alienHandler->AddRunNumber(297367);
        alienHandler->AddRunNumber(297366);
        alienHandler->AddRunNumber(297363);
        alienHandler->AddRunNumber(297336);
        alienHandler->AddRunNumber(297335);
        alienHandler->AddRunNumber(297333);
        alienHandler->AddRunNumber(297332);
        alienHandler->AddRunNumber(297317);
        alienHandler->AddRunNumber(297311);
        alienHandler->AddRunNumber(297310);
        alienHandler->AddRunNumber(297278);
        alienHandler->AddRunNumber(297222);
        alienHandler->AddRunNumber(297221);
        alienHandler->AddRunNumber(297218);
        alienHandler->AddRunNumber(297196);
        alienHandler->AddRunNumber(297195);
        alienHandler->AddRunNumber(297193);
        alienHandler->AddRunNumber(297133);
        alienHandler->AddRunNumber(297132);
        alienHandler->AddRunNumber(297129);
        alienHandler->AddRunNumber(297128);
        alienHandler->AddRunNumber(297124);
        alienHandler->AddRunNumber(297123);
        alienHandler->AddRunNumber(297119);
        alienHandler->AddRunNumber(297118);
        alienHandler->AddRunNumber(297117);
        alienHandler->AddRunNumber(297085);
        alienHandler->AddRunNumber(297035);
        alienHandler->AddRunNumber(297031);
        alienHandler->AddRunNumber(296966);
        alienHandler->AddRunNumber(296941);
        alienHandler->AddRunNumber(296938);
        alienHandler->AddRunNumber(296935);
        alienHandler->AddRunNumber(296934);

        // alienHandler->AddRunNumber(296932);
        // alienHandler->AddRunNumber(296931);
        // alienHandler->AddRunNumber(296930);
        // alienHandler->AddRunNumber(296903);
        // alienHandler->AddRunNumber(296900);
        // alienHandler->AddRunNumber(296899);
        // alienHandler->AddRunNumber(296894);
        // alienHandler->AddRunNumber(296852);
        // alienHandler->AddRunNumber(296851);
        // alienHandler->AddRunNumber(296850);
        // alienHandler->AddRunNumber(296848);
        // alienHandler->AddRunNumber(296839);
        // alienHandler->AddRunNumber(296838);
        // alienHandler->AddRunNumber(296836);
        // alienHandler->AddRunNumber(296835);
        // alienHandler->AddRunNumber(296799);
        // alienHandler->AddRunNumber(296794);
        // alienHandler->AddRunNumber(296793);
        // alienHandler->AddRunNumber(296790);
        // alienHandler->AddRunNumber(296787);
        // alienHandler->AddRunNumber(296786);
        // alienHandler->AddRunNumber(296785);
        // alienHandler->AddRunNumber(296784);
        // alienHandler->AddRunNumber(296781);
        // alienHandler->AddRunNumber(296752);
        // alienHandler->AddRunNumber(296694);
        // alienHandler->AddRunNumber(296693);
        // alienHandler->AddRunNumber(296691);
        // alienHandler->AddRunNumber(296690);
      }
    }

    // number of files per subjob
    alienHandler->SetSplitMaxInputFileNumber(40);
    alienHandler->SetExecutable("myTask.sh");
    // specify how many seconds your job may take
    alienHandler->SetTTL(30000);
    alienHandler->SetJDLName("myTask.jdl");

    alienHandler->SetOutputToRunNo(kTRUE);
    alienHandler->SetKeepLogs(kTRUE);
    // merging: run with kTRUE to merge on grid
    // after re-running the jobs in SetRunMode("terminate")
    // (see below) mode, set SetMergeViaJDL(kFALSE)
    // to collect final results
    alienHandler->SetMaxMergeStages(1);
    if (gridRun || gridMerge)
    {
      alienHandler->SetMergeViaJDL(kTRUE);
    }
    if (localMerge)
    {
      alienHandler->SetMergeViaJDL(kFALSE);
    }

    // define the output folders
    if (dataset.EqualTo("LHC10h"))
      alienHandler->SetGridWorkingDir(Form("%s/LHC10h", taskName.Data()));
    if (dataset.EqualTo("LHC15o"))
      alienHandler->SetGridWorkingDir(Form("%s/LHC15o", taskName.Data()));
    if (dataset.EqualTo("LHC18q"))
      alienHandler->SetGridWorkingDir(Form("%s/LHC18q", taskName.Data()));
    if (dataset.EqualTo("LHC18r"))
      alienHandler->SetGridWorkingDir(Form("%s/LHC18r", taskName.Data()));

    alienHandler->SetGridOutputDir("output");

    // connect the alien plugin to the manager
    mgr->SetGridHandler(alienHandler);
    if (gridTest)
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
      if (gridMerge || localMerge)
        alienHandler->SetRunMode("terminate");
      if (gridRun)
        alienHandler->SetRunMode("full");
      mgr->StartAnalysis("grid");
    }
  }
}
