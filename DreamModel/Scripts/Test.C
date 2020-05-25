#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TAxis.h"
#include "Model.h"
#include <iostream>

void Test() {
  std::cout << "Hello Buddy" << std::endl;
  const unsigned NumMomBins = 400;
  const double rMin = 0.08;
  const double rMax = 10.0;

  TString B2effplot = "B2kfir.pdf";
  TCanvas *binplot2 = new TCanvas("binplot2", "Graph Draw Options", 1000, 1000);
  binplot2->SetLeftMargin(0.1505);
  binplot2->SetRightMargin(0.035);

  Model *Gauss = new Model();

  //Gauss->SetHulthenWF(true);
  Gauss->SetMomBins(NumMomBins, rMin, rMax);
  Gauss->LoadWF();

  TGraph *grB2ChiEFT = new TGraph();
  grB2ChiEFT->SetName("grB2ChiEFT");
  grB2ChiEFT->Set(NumMomBins);
  TGraph *grB2Hulthen = new TGraph();
  grB2Hulthen->SetName("grB2Hulthen");
  grB2Hulthen->Set(NumMomBins);
  TGraph *grB2Gauss = new TGraph();
  grB2Gauss->SetName("grB2Gauss");
  grB2Gauss->Set(NumMomBins);

  TGraph *grB2Classical = new TGraph();
  grB2Classical->SetName("grB2Classical");
  grB2Classical->Set(NumMomBins);

  TLegend *leg = new TLegend(0.55, 0.54, 0.85, 0.88);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetNColumns(1);
  leg->SetLineColor(0);

  leg->AddEntry(grB2ChiEFT, "#it{B}_{2} Chiral EFT", "pef");
  leg->AddEntry(grB2Hulthen, " #it{B}_{2} Hulthen", "pef");
  //leg->AddEntry(B2mTTwoGauss_syst, "Two gaussians", "pef");
  leg->AddEntry(grB2Gauss, "#it{B}_{2} Gaussian", "pef");
  leg->AddEntry(grB2Classical, "#it{B}_{2} Classical", "pef");

  for (unsigned uMom = 0; uMom < NumMomBins; uMom++) {
    grB2ChiEFT->SetPoint(uMom, Gauss->GetMomentum(uMom),
                         Gauss->GetB2ChiEFT(uMom));  //
    grB2Hulthen->SetPoint(uMom, Gauss->GetMomentum(uMom),
                          Gauss->GetB2Hulthen(uMom));  //
    grB2Gauss->SetPoint(uMom, Gauss->GetMomentum(uMom),
                        Gauss->GetB2_kfir(uMom));  //
    grB2Classical->SetPoint(uMom, Gauss->GetMomentum(uMom),
                        Gauss->GetB2_Classical(uMom));  //
    std::cout<<"B2_hulthen = "<< Gauss->GetB2Hulthen(uMom)<< "Radius ="<< Gauss->GetMomentum(uMom)<< std::endl;
  }
  grB2Hulthen->SetTitle("; #it{R} (fm); #it{B}_{2} (GeV^{2}/c^{3})");
  grB2Hulthen->GetXaxis()->SetLimits(0.08, 10.1);//our concern is R from 0.8 to 5 fm as we know that this range is most useful to study coalescence (available source size falls within this range :P)
  grB2Hulthen->GetYaxis()->SetRangeUser(1e-5, 1e3);

  grB2Hulthen->SetMarkerColor(65 + 2);
  grB2Hulthen->SetLineColor(65 + 2);

  grB2Gauss->SetMarkerColor(kYellow + 2);
  grB2Gauss->SetLineColor(kYellow + 2);

  grB2ChiEFT->SetMarkerColor(kBlue + 2);
  grB2ChiEFT->SetLineColor(kBlue + 2);

  grB2Classical->SetMarkerColor(kBlack + 0);
  grB2Classical->SetLineColor(kBlack + 0);

  grB2Hulthen->Draw("ALP");
  grB2Gauss->Draw("same");
  grB2ChiEFT->Draw("same");
  grB2Classical->Draw("same");
  leg->Draw("same");

  gPad->SetLogy();
  gPad->SetLogx();

  binplot2->SaveAs(B2effplot);
  binplot2->Close();
  delete Gauss;
  TString nameWFplot = "WFplot.pdf";
  TCanvas *WFplot = new TCanvas("WFplot", "Graph Draw Options", 1000, 1000);
  WFplot->SetLeftMargin(0.1505);
  WFplot->SetRightMargin(0.035);
  const unsigned NumrBins = 400;
  const double r_Min = 0.0;
  const double r_Max = 10;


  Model *M1 = new Model();
  Model *M2 = new Model();
  Model *M3 = new Model();

  M1->SetMomBins(NumrBins, r_Min, r_Max);
  M1->LoadWF();
  M2->SetMomBins(NumrBins, r_Min, r_Max);
  M2->LoadWF();
  M3->SetMomBins(NumrBins, r_Min, r_Max);
  M3->LoadWF();
  M1->SetGaussWF(true);
  M2->SetHulthenWF(true);
  M3->SetChiEFTWF(true);
  M1->SetDrawWF(true);
  M2->SetDrawWF(true);
  M3->SetDrawWF(true);

  TGraph *grWFChiEFT = new TGraph();
  grWFChiEFT->SetName("grWFChiEFT");
  grWFChiEFT->Set(NumrBins);
  TGraph *grWFHulthen = new TGraph();
  grWFHulthen->SetName("grWFHulthen");
  grWFHulthen->Set(NumrBins);
  TGraph *grWFGauss = new TGraph();
  grWFGauss->SetName("grWFGauss");
  grWFGauss->Set(NumrBins);
  TLegend *leg2= new TLegend(0.5, 0.6, 0.8, 0.87);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.040);
    leg2->SetLineColor(0);
    leg2->AddEntry(grWFGauss, "Gauss", "pef");
    leg2->AddEntry(grWFHulthen, "Hulthen", "pef");
    leg2->AddEntry(grWFChiEFT, "Chiral EFT", "pef");
  for (unsigned uMom = 0; uMom < NumrBins; uMom++) {
    grWFChiEFT->SetPoint(uMom, M3->GetMomentum(uMom),
                         M3->GetWaveFunction(uMom));  //
    grWFHulthen->SetPoint(uMom, M2->GetMomentum(uMom),
                          M2->GetWaveFunction(uMom));  //
    grWFGauss->SetPoint(uMom, M1->GetMomentum(uMom), M1->GetWaveFunction(uMom));  //
  }
  grWFGauss->SetTitle("; #it{r} (fm); #psi");
  grWFGauss->GetXaxis()->SetLimits(0, 11);
  grWFGauss->GetYaxis()->SetRangeUser(-0.01, 0.35);

  grWFGauss->SetMarkerColor(kYellow + 2);
  grWFGauss->SetLineColor(kYellow + 2);

  grWFHulthen->SetMarkerColor(65 + 2);
  grWFHulthen->SetLineColor(65 + 2);

  grWFChiEFT->SetMarkerColor(kBlue + 2);
  grWFChiEFT->SetLineColor(kBlue + 2);

  grWFGauss->Draw("ALP");
  grWFHulthen->Draw("same");
  grWFChiEFT->Draw("same");
  leg2->Draw("same");
  //gPad->SetLogy();
  //gPad->SetLogx();
  WFplot->SaveAs(nameWFplot);
  WFplot->Close();
  delete M1;
  delete M2;
  delete M3;

}
