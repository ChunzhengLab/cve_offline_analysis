TChain* CreateChain(const char *xmlfile, const char *type="ESD");

const char *anatype = "AOD";

void lego_train()
{
// Analysis using AOD data
// Automatically generated analysis steering macro executed in grid subjobs

   TStopwatch timer;
   timer.Start();

// Set temporary merging directory to current one
   gSystem->Setenv("TMPDIR", gSystem->pwd());

// Set temporary compilation directory to current one
   gSystem->SetBuildDir(gSystem->pwd(), kTRUE);

// Reset existing include path and add current directory first in the search
   gSystem->SetIncludePath("-I.");

// analysis source to be compiled at runtime (if any)

// read the analysis manager from file
   AliAnalysisManager *mgr = AliAnalysisAlien::LoadAnalysisManager("lego_train.root");
   if (!mgr) return;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   plugin->SetRunMode("test");
   plugin->SetFileForTestMode("data.txt");
   plugin->SetNtestFiles(2);
   mgr->SetGridHandler(plugin);
   mgr->SetDebugLevel(0);
   mgr->SetNSysInfo(40);
   mgr->PrintStatus();
   AliLog::SetGlobalLogLevel(AliLog::kError);
   mgr->StartAnalysis("localfile", 123456789, 0);
   timer.Stop();
   timer.Print();
}

