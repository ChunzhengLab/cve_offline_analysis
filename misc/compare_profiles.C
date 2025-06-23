void compare_profiles() {
    // Helper function to get profiles from a directory
    auto getProfiles = [](TFile* file, const char* dirName, const char* listName,
                         const char* profileNameBase = "DeltaLambdaProtonMass") {
        struct ProfilePair {
            TProfile3D* prof_0 = nullptr;
            TProfile3D* prof_1 = nullptr;
        };

        ProfilePair result;

        TDirectory* dir = file->GetDirectory(dirName);
        if (!dir) {
            std::cerr << "Error: " << dirName << " directory not found in file" << std::endl;
            return result;
        }

        TList* list = nullptr;
        dir->GetObject(listName, list);
        if (!list) {
            std::cerr << "Error: " << listName << " not found in directory" << std::endl;
            return result;
        }

        // Construct profile names based on the base name
        TString diffName_0 = Form("fProfile3DDiff%sDEta_0", profileNameBase);
        TString name_0 = Form("fProfile3D%sDEta_0", profileNameBase);
        TString diffName_1 = Form("fProfile3DDiff%sDEta_1", profileNameBase);
        TString name_1 = Form("fProfile3D%sDEta_1", profileNameBase);

        // Print out what we're looking for if it's Gamma
        if (strcmp(profileNameBase, "GammaLambdaProtonMass") == 0) {
            std::cout << "Searching for Gamma profiles: " << diffName_0.Data() << " or " << name_0.Data() << std::endl;
        }

        // Try both possible profile names
        result.prof_0 = dynamic_cast<TProfile3D*>(list->FindObject(diffName_0));
        if (!result.prof_0) {
            result.prof_0 = dynamic_cast<TProfile3D*>(list->FindObject(name_0));
        }

        result.prof_1 = dynamic_cast<TProfile3D*>(list->FindObject(diffName_1));
        if (!result.prof_1) {
            result.prof_1 = dynamic_cast<TProfile3D*>(list->FindObject(name_1));
        }

        // Print result for Gamma profiles
        if (strcmp(profileNameBase, "GammaLambdaProtonMass") == 0) {
            if (result.prof_0) {
                std::cout << "Found Gamma profile _0 in " << dirName << std::endl;
            } else {
                std::cout << "WARNING: Gamma profile _0 not found in " << dirName << std::endl;
            }

            if (result.prof_1) {
                std::cout << "Found Gamma profile _1 in " << dirName << std::endl;
            } else {
                std::cout << "WARNING: Gamma profile _1 not found in " << dirName << std::endl;
            }
        }

        return result;
    };

    // Convert 3D profile to 1D profile
    auto convertTo1D = [](TProfile3D* prof3D, const char* name) -> TProfile* {
        if (!prof3D) return nullptr;

        TProfile2D* prof2D = prof3D->Project3DProfile("xz");
        return prof2D->ProfileY(name, 15, 16);
    };

    // Draw all profiles together
    auto drawAllProfiles = [](TVirtualPad* pad,
                            TProfile* default_0_old, TProfile* default_1_old,
                            TProfile* nue_0_old, TProfile* nue_1_old,
                            TProfile* default_0_new, TProfile* default_1_new) {
        pad->cd();
        pad->DrawFrame(0, -0.01, 60, 0.01);

        TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
        leg->SetHeader("Profiles Comparison");

        bool firstDrawn = false;

        // Draw default_old profiles (blue)
        if (default_0_old) {
            default_0_old->SetLineColor(kBlue);
            default_0_old->SetLineWidth(2);
            default_0_old->SetTitle("TProfile X at ybin=15,16;X axis;Profile value");
            default_0_old->Draw("EL same");
            leg->AddEntry(default_0_old, "default_old_frame", "l");
            firstDrawn = true;
        }

        if (default_1_old) {
            default_1_old->SetLineColor(kBlue);
            default_1_old->SetLineWidth(2);
            default_1_old->SetLineStyle(2); // Dashed line for _1 profiles

            if (firstDrawn) {
                default_1_old->Draw("L SAME");
            } else {
                default_1_old->SetTitle("TProfile X at ybin=15,16;X axis;Profile value");
                default_1_old->Draw("EL same");
                firstDrawn = true;
            }
        }

        // Draw NUE profiles (red)
        if (nue_0_old) {
            nue_0_old->SetLineColor(kRed);
            nue_0_old->SetLineWidth(2);

            if (firstDrawn) {
                nue_0_old->Draw("L SAME");
            } else {
                nue_0_old->SetTitle("TProfile X at ybin=15,16;X axis;Profile value");
                nue_0_old->Draw("EL same");
                firstDrawn = true;
            }

            leg->AddEntry(nue_0_old, "doLambdaNUE_old_frame", "l");
        }

        if (nue_1_old) {
            nue_1_old->SetLineColor(kRed);
            nue_1_old->SetLineWidth(2);
            nue_1_old->SetLineStyle(2); // Dashed line for _1 profiles

            if (firstDrawn) {
                nue_1_old->Draw("L SAME");
            } else {
                nue_1_old->SetTitle("TProfile X at ybin=15,16;X axis;Profile value");
                nue_1_old->Draw("EL same");
                firstDrawn = true;
            }
        }

        // Draw new_frame profiles (black)
        if (default_0_new) {
            default_0_new->SetLineColor(kBlack);
            default_0_new->SetLineWidth(2);

            if (firstDrawn) {
                default_0_new->Draw("L SAME");
            } else {
                default_0_new->SetTitle("TProfile X at ybin=15,16;X axis;Profile value");
                default_0_new->Draw("EL same");
                firstDrawn = true;
            }

            leg->AddEntry(default_0_new, "new_frame", "l");
        }

        if (default_1_new) {
            default_1_new->SetLineColor(kBlack);
            default_1_new->SetLineWidth(2);
            default_1_new->SetLineStyle(2); // Dashed line for _1 profiles

            if (firstDrawn) {
                default_1_new->Draw("L SAME");
            } else {
                default_1_new->SetTitle("TProfile X at ybin=15,16;X axis;Profile value");
                default_1_new->Draw("EL same");
            }
        }

        leg->Draw();
    };

    // Open files
    TFile* file_old = TFile::Open("AnalysisResults_old.root");
    if (!file_old || file_old->IsZombie()) {
        std::cerr << "Error: Cannot open AnalysisResults_old.root" << std::endl;
        return;
    }

    TFile* file_new = TFile::Open("AnalysisResults_new.root");
    if (!file_new || file_new->IsZombie()) {
        std::cerr << "Error: Cannot open AnalysisResults_new.root" << std::endl;
        return;
    }

    // Get 3D profiles from old file
    auto profiles_default_old = getProfiles(file_old, "default", "ListResults_default");
    auto profiles_nue_old = getProfiles(file_old, "DoLambdaNUE", "ListResults_DoLambdaNUE");

    // Get 3D profiles from new file (only default directory)
    auto profiles_default_new = getProfiles(file_new, "default", "ListResults_default");

    // Check if at least one profile is available
    if (!profiles_default_old.prof_0 && !profiles_default_old.prof_1 &&
        !profiles_nue_old.prof_0 && !profiles_nue_old.prof_1 &&
        !profiles_default_new.prof_0 && !profiles_default_new.prof_1) {
        std::cerr << "Error: No profiles available to plot" << std::endl;
        return;
    }

    // Convert 3D profiles to 1D profiles for plotting
    TProfile* profile_default_0_old = convertTo1D(profiles_default_old.prof_0, "profile_default_0_old");
    TProfile* profile_default_1_old = convertTo1D(profiles_default_old.prof_1, "profile_default_1_old");
    TProfile* profile_nue_0_old = convertTo1D(profiles_nue_old.prof_0, "profile_nue_0_old");
    TProfile* profile_nue_1_old = convertTo1D(profiles_nue_old.prof_1, "profile_nue_1_old");

    TProfile* profile_default_0_new = convertTo1D(profiles_default_new.prof_0, "profile_default_0_new");
    TProfile* profile_default_1_new = convertTo1D(profiles_default_new.prof_1, "profile_default_1_new");

    // Get Gamma profiles from old file
    auto gamma_profiles_default_old = getProfiles(file_old, "default", "ListResults_default", "GammaLambdaProtonMass");
    auto gamma_profiles_nue_old = getProfiles(file_old, "DoLambdaNUE", "ListResults_DoLambdaNUE", "GammaLambdaProtonMass");

    // Get Gamma profiles from new file
    auto gamma_profiles_default_new = getProfiles(file_new, "default", "ListResults_default", "GammaLambdaProtonMass");

    // Convert 3D profiles to 1D profiles for plotting (Gamma)
    TProfile* gamma_profile_default_0_old = convertTo1D(gamma_profiles_default_old.prof_0, "gamma_profile_default_0_old");
    TProfile* gamma_profile_default_1_old = convertTo1D(gamma_profiles_default_old.prof_1, "gamma_profile_default_1_old");
    TProfile* gamma_profile_nue_0_old = convertTo1D(gamma_profiles_nue_old.prof_0, "gamma_profile_nue_0_old");
    TProfile* gamma_profile_nue_1_old = convertTo1D(gamma_profiles_nue_old.prof_1, "gamma_profile_nue_1_old");

    TProfile* gamma_profile_default_0_new = convertTo1D(gamma_profiles_default_new.prof_0, "gamma_profile_default_0_new");
    TProfile* gamma_profile_default_1_new = convertTo1D(gamma_profiles_default_new.prof_1, "gamma_profile_default_1_new");

    // Create canvas and divide into two parts
    TCanvas* c = new TCanvas("c", "Compare profiles", 1800, 800);
    c->Divide(2, 1);

    // Draw Lambda profiles on the left pad
    TVirtualPad* pad1 = c->cd(1);
    drawAllProfiles(pad1,
                  profile_default_0_old, profile_default_1_old,
                  profile_nue_0_old, profile_nue_1_old,
                  profile_default_0_new, profile_default_1_new);

    // Draw Gamma profiles on the right pad
    TVirtualPad* pad2 = c->cd(2);
    drawAllProfiles(pad2,
                  gamma_profile_default_0_old, gamma_profile_default_1_old,
                  gamma_profile_nue_0_old, gamma_profile_nue_1_old,
                  gamma_profile_default_0_new, gamma_profile_default_1_new);

    c->Update();
}
