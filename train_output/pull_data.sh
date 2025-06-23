#!/bin/bash

# This script pulls AnalysisResults.root files from ALIEN storage for different ALICE runs and particles.
# It uses alien_cp to copy the files to the local directory.

set -e
alien_cp alien::/alice/cern.ch/user/a/alitrain/PWGCF/CF_PbPb/16716_20250601-1912/merge/AnalysisResults.root file:./AnalysisResults_CVE2025_18q_TPC_Proton.root;
sleep 2
echo "18q_TPC_Proton.root downloaded" #

alien_cp alien::/alice/cern.ch/user/a/alitrain/PWGCF/CF_PbPb/16718_20250601-1912/merge/AnalysisResults.root file:./AnalysisResults_CVE2025_18q_TPC_Pion.root;
sleep 2
echo "18q_TPC_Pion.root downloaded" #

alien_cp alien::/alice/cern.ch/user/a/alitrain/PWGCF/CF_PbPb/16717_20250601-1912/merge/AnalysisResults.root file:./AnalysisResults_CVE2025_18q_TPC_Lambda.root;
sleep 2
echo "18q_TPC_Lambda.root downloaded" #

alien_cp alien::/alice/cern.ch/user/a/alitrain/PWGCF/CF_PbPb/16715_20250601-1912/merge/AnalysisResults.root file:./AnalysisResults_CVE2025_18q_TPC_Hadron.root;
sleep 2
echo "18q_TPC_Hadron.root downloaded" #

alien_cp alien::/alice/cern.ch/user/a/alitrain/PWGCF/CF_PbPb/16722_20250601-1912/merge/AnalysisResults.root file:./AnalysisResults_CVE2025_18r_TPC_Proton.root;
sleep 2
echo "18r_TPC_Proton.root downloaded" #

alien_cp alien::/alice/cern.ch/user/a/alitrain/PWGCF/CF_PbPb/16720_20250601-1912/merge/AnalysisResults.root file:./AnalysisResults_CVE2025_18r_TPC_Lambda.root;
sleep 2
echo "18r_TPC_Lambda.root downloaded" #75.8%

alien_cp alien::/alice/cern.ch/user/a/alitrain/PWGCF/CF_PbPb/16719_20250601-1912/merge/AnalysisResults.root file:./AnalysisResults_CVE2025_18r_TPC_Hadron.root;
sleep 2
echo "18r_TPC_Hadron.root downloaded" #77.4%

alien_cp alien::/alice/cern.ch/user/a/alitrain/PWGCF/CF_PbPb/16721_20250601-1912/merge/AnalysisResults.root file:./AnalysisResults_CVE2025_18r_TPC_Pion.root;
sleep 2
echo "18r_TPC_Pion.root downloaded" #85.1%
