#ifndef ALI_TRAIN_GENERATOR_H
#define ALI_TRAIN_GENERATOR_H

#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCfg.h"
#include "TGrid.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

// —— 辅助：读取环境变量 ——
static Int_t GetEnvInt(const char *name, Int_t def = 0) {
  const char *v = gSystem->Getenv(name);
  return v ? TString(v).Atoi() : def;
}
static TString GetEnvStr(const char *name, const char *def = "") {
  const char *v = gSystem->Getenv(name);
  return v ? TString(v) : TString(def);
}

// —— 辅助：收集 TEST_DIR 及其子目录 ——
static std::vector<TString> CollectTestDirs() {
  std::vector<TString> dirs;
  dirs.push_back(GetEnvStr("TEST_DIR"));
  for (int i = 1;; ++i) {
    TString child = GetEnvStr(Form("TEST_DIR_child_%d", i), "-1");
    if (child == "-1")
      break;
    dirs.push_back(child);
  }
  return dirs;
}

// —— 辅助：配置 AliAnalysisAlien 插件 ——
static void ConfigurePlugin(AliAnalysisAlien *plugin) {
  plugin->SetProductionMode();
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAdditionalRootLibs(
      "libVMC.so libPhysics.so libTree.so libMinuit.so "
      "libProof.so libSTEERBase.so libESD.so libAOD.so");
  plugin->SetJobTag("test/test");
  plugin->SetMaxMergeFiles(GetEnvInt("MAX_MERGE_FILES", 10));
  plugin->SetTTL(GetEnvInt("TTL", 3600));
  plugin->SetAnalysisMacro("lego_train.C");
  plugin->SetValidationScript("validation.sh");

  // 排除文件列表
  TString excl = GetEnvStr("EXCLUDE_FILES");
  plugin->SetRegisterExcludes(excl + " AliAOD.root");

  // —— 友链配置 ——
  TString names = GetEnvStr("FRIEND_CHAIN_NAMES");
  TString libs = GetEnvStr("FRIEND_CHAIN_LIBRARIES");
  if (!names.IsNull()) {
    if (!libs.IsNull()) {
      plugin->SetFriendChainName(names, libs);
    } else {
      plugin->SetFriendChainName(names);
    }
  } else if (!libs.IsNull()) {
    // 先取出已有的 rootLibs，拼接后再设置
    TString allLibs = plugin->GetAdditionalRootLibs();
    allLibs += " ";
    allLibs += libs;
    plugin->SetAdditionalRootLibs(allLibs);
  }

  // 外部包
  TString extra = GetEnvStr("ADDITIONAL_PACKAGES");
  plugin->AddExternalPackage(extra + " jemalloc::v3.6.0");

  // JDL & 脚本
  plugin->SetJDLName("lego_train.jdl");
  plugin->SetExecutable("lego_train.sh");
  plugin->SetSplitMode("se");
  plugin->SetGridOutputDir("./");
  plugin->SetGridWorkingDir("./");

  // 日志保留由调用处决定
  plugin->SetMergeViaJDL();
}

// —— 辅助：添加数据文件 ——
static void AddDataFiles(AliAnalysisAlien *plugin,
                         const std::vector<TString> &dirs,
                         const TString &anchor, Int_t nFiles, Int_t AOD,
                         const TString &periodName) {
  TString archive = (AOD == 2 ? "aod_archive.zip" : "root_archive.zip");
  TString friendNames = GetEnvStr("FRIEND_CHAIN_NAMES");
  if (!friendNames.IsNull()) {
    friendNames.ReplaceAll(" ", ";");
    archive += ";" + friendNames;
  }
  plugin->SetNtestFiles(dirs.size() * nFiles);

  if (AOD == 3) {
    Printf(">>>> Expecting MC only production");
    plugin->SetUseMCchain();
  }

  for (auto &dir : dirs) {
    Bool_t special =
        (periodName == "AMPT_LHC12g6" || periodName == "AMPT_LHC12c3");
    if (special) {
      plugin->AddDataFile(Form("%s/Kinematics.root", dir.Data()));
      plugin->AddDataFile(Form("%s/%s", dir.Data(), anchor.Data()));
    } else {
      TString path = Form("%s/%s", dir.Data(), anchor.Data());
      Printf("Adding dataset path: %s", path.Data());
      plugin->AddDataFile(path);
    }
  }
}

// —— 辅助：设置执行命令与分割策略 ——
static void SetupExecution(AliAnalysisAlien *plugin, Int_t AOD, Int_t splitMax,
                           Int_t nTestEvents) {
  if (AOD == 100) {
    Long64_t totalEvents = TString(gSystem->Getenv("GEN_TOTAL_EVENTS")).Atoll();
    Long64_t jobs = totalEvents / splitMax;
    plugin->SetMCLoop(true);
    plugin->SetSplitMode(Form("production:1-%lld", jobs));
    plugin->SetNMCjobs(jobs);
    plugin->SetNMCevents(nTestEvents);
    plugin->SetExecutableCommand("aliroot -b -q");
  } else {
    plugin->SetSplitMaxInputFileNumber(splitMax);
    plugin->SetExecutableCommand("root -b -q");
    plugin->SetInputFormat("xml-single");
  }
}

// —— 主入口：生成 TRAIN 或 TEST ——
void generate(const char *module = "__ALL__") {
  // 环境变量读取
  Int_t nFiles = 2;
  Int_t nTestEvents = GetEnvInt("TEST_FILES_NO", 2);
  TString anchor = GetEnvStr("FILE_PATTERN", "AliAOD.root");
  Int_t splitMax = GetEnvInt("SPLIT_MAX_INPUT_FILE_NUMBER", 1);
  Int_t maxMerge = GetEnvInt("MAX_MERGE_FILES", 10);
  Int_t debugLevel = GetEnvInt("DEBUG_LEVEL", 0);
  Int_t ttl = GetEnvInt("TTL", 3600);
  Int_t AOD = GetEnvInt("AOD", 1);
  Bool_t isPP = (strcmp(gSystem->Getenv("PP"), "true") == 0);
  TString periodName = GetEnvStr("PERIOD_NAME", "");
  TString excl = GetEnvStr("EXCLUDE_FILES", "");
  TString outFiles = GetEnvStr("OUTPUT_FILES", "");

  // 模式解析
  Bool_t doProd = kFALSE;
  if (strcmp(module, "__TRAIN__") == 0) {
    doProd = kTRUE;
    module = "";
  } else if (strcmp(module, "__ALL__") == 0) {
    module = "";
  }

  // AliEn 连接
  if (GetEnvInt("ADDTASK_NEEDS_ALIEN", 0) == 1) {
    Printf("Connecting to AliEn...");
    TGrid::Connect("alien:");
  }

  // 加载任务配置
  TObjArray *tasks =
      AliAnalysisTaskCfg::ExtractModulesFrom("MLTrainDefinition.cfg");
  Printf(">>>>>>> Read train configuration");
  tasks->Print();

  // 创建并配置 plugin
  AliAnalysisAlien *plugin = new AliAnalysisAlien("lego_train");
  ConfigurePlugin(plugin);
  plugin->SetMaxMergeFiles(maxMerge);
  plugin->SetTTL(ttl);
  if (!doProd)
    plugin->SetKeepLogs(kTRUE);

  // 添加数据
  std::vector<TString> testDirs = CollectTestDirs();
  AddDataFiles(plugin, testDirs, anchor, nFiles, AOD, periodName);

  // 全局变量宏
  Int_t err = 0;
  gROOT->Macro("globalvariables.C", &err);
  if (err != 0) {
    Printf("ERROR: globalvariables.C failed. Exiting.");
    return;
  }

  // 加载模块并创建管理器
  plugin->AddModules(tasks);
  plugin->CreateAnalysisManager("train", "handlers.C");

  // 执行命令与切分
  SetupExecution(plugin, AOD, splitMax, nTestEvents);

  // 调整 AnalysisManager
  auto mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->SetDebugLevel(debugLevel);
  mgr->SetNSysInfo(isPP ? 1000 : 40);
  mgr->SetFileInfoLog("fileinfo.log");

  // 生成
  if (doProd)
    plugin->GenerateTrain("lego_train");
  else
    plugin->GenerateTest("lego_train", module);

  // 输出合法性校验
  TString valid = outFiles + "," + excl;
  TString actual = plugin->GetListOfFiles("out");
  TObjArray *toks = actual.Tokenize(",");
  Bool_t ok = kTRUE;
  for (Int_t i = 0; i < toks->GetEntries(); ++i) {
    TString f = toks->At(i)->GetName();
    if (!valid.Contains(f)) {
      Printf("ERROR: Output file %s not allowed (should be in %s)", f.Data(),
             valid.Data());
      ok = kFALSE;
    }
  }
  delete toks;
  if (!ok) {
    Printf(">>>> Invalid output files requested. Deleting lego_train.C");
    gSystem->Unlink("lego_train.C");
  }
}

#endif // ALI_TRAIN_GENERATOR_H
