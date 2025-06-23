void plot_QVector() {
    auto file = TFile::Open("/Users/wangchunzheng/works/Experiments/cve/train_output/AnalysisResults_CVE2025_18q_TPC_Hadron.root");
    auto list = (TList*)file->Get("default/ListQA_default");
    auto p2_Qx_bfRC = (TProfile2D*)list->FindObject("fProfile2DQxCentVz_bfRC");
    auto p2_Qy_bfRC = (TProfile2D*)list->FindObject("fProfile2DQyCentVz_bfRC");
    auto p_Qx_bfRC = (TProfile*)p2_Qx_bfRC->ProjectionX("h_Qx_bfRC", 1, 3);
    auto p_Qy_bfRC = (TProfile*)p2_Qy_bfRC->ProjectionX("h_Qy_bfRC", 1, 3);

    auto p2_Qx_afRC = (TProfile2D*)list->FindObject("fProfile2DQxCentVz_afRC");
    auto p2_Qy_afRC = (TProfile2D*)list->FindObject("fProfile2DQyCentVz_afRC");
    auto p_Qx_afRC = (TProfile*)p2_Qx_afRC->ProjectionX("h_Qx_afRC", 1, 3);
    auto p_Qy_afRC = (TProfile*)p2_Qy_afRC->ProjectionX("h_Qy_afRC", 1, 3);

    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->Divide(2, 1);

    canvas->cd(1)->DrawFrame(0.,-0.008,70.,0.008,"TPC Q-vector Before Recenter;Centrality;Q-vector");
    p_Qy_bfRC->SetLineColor(kRed);
    p_Qy_bfRC->SetMarkerStyle(kFullCircle);
    p_Qy_bfRC->SetMarkerSize(1);
    p_Qy_bfRC->SetMarkerColor(kRed);
    p_Qy_bfRC->Draw("same");
    p_Qx_bfRC->SetLineColor(kBlue);
    p_Qx_bfRC->SetMarkerStyle(kFullSquare);
    p_Qx_bfRC->SetMarkerSize(1);
    p_Qx_bfRC->SetMarkerColor(kBlue);
    p_Qx_bfRC->Draw("same");
    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(p_Qy_bfRC, "Qy", "lp");
    legend->AddEntry(p_Qx_bfRC, "Qx", "lp");
    legend->Draw();

    canvas->cd(2)->DrawFrame(0.,-0.008,70.,0.008,"TPC Q-vector After Recenter;Centrality;Q-vector");
    p_Qy_afRC->SetLineColor(kRed);
    p_Qy_afRC->SetMarkerStyle(kFullCircle);
    p_Qy_afRC->SetMarkerSize(1);
    p_Qy_afRC->SetMarkerColor(kRed);
    p_Qy_afRC->Draw("same");
    p_Qx_afRC->SetLineColor(kBlue);
    p_Qx_afRC->SetMarkerStyle(kFullSquare);
    p_Qx_afRC->SetMarkerSize(1);
    p_Qx_afRC->SetMarkerColor(kBlue);
    p_Qx_afRC->Draw("same");
    legend->Draw();
}
