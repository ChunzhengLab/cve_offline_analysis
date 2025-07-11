void lego_train_merge(const char *dir, Int_t stage=0)
{
// Automatically generated merging macro executed in grid subjobs

   TStopwatch timer;
   timer.Start();

// Reset existing include path and add current directory first in the search
   gSystem->SetIncludePath("-I.");

// Analysis source to be compiled at runtime (if any)

// Connect to AliEn
   if (!TGrid::Connect("alien://")) return;
// Set temporary merging directory to current one
   gSystem->Setenv("TMPDIR", gSystem->pwd());

// Set temporary compilation directory to current one
   gSystem->SetBuildDir(gSystem->pwd(), kTRUE);

   TString outputDir = dir;
   TString outputFiles = "AnalysisResults.root,EventStat_temp.root";
   TString mergeExcludes = " QnEventQA.root CalibrationQA.root QnAnalysisEventQA.root AliAOD.root";
   TObjArray *list = outputFiles.Tokenize(",");
   TIter *iter = new TIter(list);
   TObjString *str;
   TString outputFile;
   Bool_t merged = kTRUE;
   AliAnalysisManager *mgr = AliAnalysisAlien::LoadAnalysisManager("lego_train.root");
   if (!mgr) {
      printf("ERROR: Analysis manager could not be extracted from file ");
      return;
   }
   printf("===================================\n");
   printf("Testing merging...\n");
   printf("===================================\n");
   while((str=(TObjString*)iter->Next())) {
      outputFile = str->GetString();
      if (outputFile.Contains("*")) continue;
      Int_t index = outputFile.Index("@");
      if (index > 0) outputFile.Remove(index);
      // Skip already merged outputs
      if (!gSystem->AccessPathName(outputFile)) {
         printf("Output file <%s> found. Not merging again.\n",outputFile.Data());
         continue;
      }
      if (mergeExcludes.Contains(outputFile.Data())) continue;
      merged = AliAnalysisAlien::MergeOutput(outputFile, outputDir, 10, stage);
      if (!merged) {
         printf("ERROR: Cannot merge %s\n", outputFile.Data());
         return;
      }
   }
   // all outputs merged, validate
   ofstream out;
   out.open("outputs_valid", ios::out);
   out.close();
   // read the analysis manager from file
   printf("===================================\n");
   printf("Testing Terminate()...\n");
   printf("===================================\n");
   mgr->SetRunFromPath(mgr->GetRunFromAlienPath(dir));
   mgr->SetSkipTerminate(kFALSE);
   mgr->PrintStatus();
   AliLog::SetGlobalLogLevel(AliLog::kError);
   TTree *tree = NULL;
   mgr->StartAnalysis("gridterminate", tree);
}

