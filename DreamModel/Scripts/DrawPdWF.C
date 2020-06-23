#include "DLM_Histo.h"
#include "DLM_CppTools.h"
#include "DLM_CppTools.cpp"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "Model.h"

#include <iostream>

void ReadPDeuteronWF(const char *InputFileName, DLM_Histo<float> &OutputU,
                     DLM_Histo<float> &OutputW) {

  FILE *InFile;
  InFile= fopen(InputFileName, "r");
  if (!InFile) {
    printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n",
           InputFileName);
    return;
  }

  char *cdummy= new char[1024];
  unsigned CurrentLine= 0;
  float fRadius, fReWF, fImWF;
  double *Radius= new double[1024];
  double *RadiusBins= new double[1025];

  unsigned NumBins= 0;
  RadiusBins[0]= 0;
  //read line-by-line until the end of the file
  while (!feof(InFile)) {
    CurrentLine++;
    //read the next line
    if (!fgets(cdummy, 1023, InFile))
      continue;
    if (CurrentLine <= 17)
      continue;
    //printf("Line#%u: %s\n",CurrentLine,cdummy);
    sscanf(cdummy, "%f %f %f", &fRadius, &fReWF, &fImWF);
    //printf("Line#%u: %.2e; %.2e; %.2e;\n",CurrentLine,fRadius,fReWF,fImWF);

    //if((CurrentLine)%int(NumRadBins)==0){
    //   sscanf(cdummy, "%f",&fMomentum);
    Radius[NumBins]= fRadius;
    if (NumBins) {
      //set the bin range in between the last two bin centers
      RadiusBins[NumBins]= 0.5 * (Radius[NumBins] + Radius[NumBins - 1]);
    }
    //}

    NumBins++;
  }
  fclose(InFile);

  //set the upper edge of the last bin, where we just add the bin width of the last bin
  //i.e. if we have l(low) c(center) u(up), we have that u=c+(c-l)=2c-l
  RadiusBins[NumBins]= 2. * Radius[NumBins - 1] - RadiusBins[NumBins - 1];

  //for(unsigned uBin=0; uBin<NumBins; uBin++){
  //    printf("#%u: %f %f%f\n",uBin,RadiusBins[uBin],Radius[uBin],RadiusBins[uBin+1]);
  //}

  //dimension (1D)
  OutputU.SetUp(1);
  OutputU.SetUp(0, NumBins, RadiusBins, Radius);
  OutputU.Initialize();

  OutputW.SetUp(1);
  OutputW.SetUp(0, NumBins, RadiusBins, Radius);
  OutputW.Initialize();

  CurrentLine= 0;
  NumBins= 0;
  InFile= fopen(InputFileName, "r");
  while (!feof(InFile)) {
    CurrentLine++;
    //read the next line
    if (!fgets(cdummy, 1023, InFile))
      continue;
    if (CurrentLine <= 17)
      continue;
    sscanf(cdummy, "%f %f %f", &fRadius, &fReWF, &fImWF);
    OutputU.SetBinContent(NumBins, fReWF);
    OutputW.SetBinContent(NumBins, fImWF);
    NumBins++;
  }
  fclose(InFile);

  delete[] cdummy;
  delete[] Radius;
}

DLM_Histo<float> hReWF;
DLM_Histo<float> hImWF;
float Evaluate_d_u(double Radius) {
  return hReWF.Eval(&Radius);
}
float Evaluate_d_w(double Radius) {
  return hImWF.Eval(&Radius);
}

void TestDLMhistoPD() {

  ReadPDeuteronWF(
      "/home/sbhawani/Desktop/PdCFTheory/WaveFunctions/Wf-pd-D/Wf-pd-D-k100-L400-l1p0.dat",
      hReWF, hImWF);

  const double RAD= 5.05;
  printf("u(%.2f) = %.4f\n", RAD, Evaluate_d_u(RAD));
  printf(" w(%.2f) = %.4f\n", RAD, Evaluate_d_w(RAD));

}

//////////
Double_t WaveChiEFT(Double_t r) {

  double u= Evaluate_d_u(r);
  double w= Evaluate_d_w(r);
  return 1 / (pow(4 * TMath::Pi() * r * r, 0.5)) * u;

}

void DrawPdWF() {
  TestDLMhistoPD();

  std::cout << "Hello Buddy" << std::endl;
  Double_t w = 900;
   Double_t h = 600;
  TString nameWFplot= "WFpdplotk100L400l1p0.pdf";
  TCanvas *WFplot= new TCanvas("WFplot", "WFplot", w, h);
  WFplot->SetLeftMargin(0.1505);
  WFplot->SetRightMargin(0.035);
  WFplot->Divide(2, 2, 0, 0);
  TPaveLabel* title = new TPaveLabel(0.1,0.94,0.9,0.97,"k = 100MeV, #Lambda = 400, Wave l1p0");
  title->Draw();
  const unsigned NumrBins= 101;
  const double r_Min= 0.5;
  const double r_Max= 97.0;

  Model *M1= new Model();
  // Model *M2 = new Model();
  // Model *M3 = new Model();

  M1->SetMomBins(NumrBins, r_Min, r_Max);
  // M1->LoadWF();
  //M1->SetGaussWF(true);
  //M1->SetDrawWF(true);

  TGraph *grRelWF= new TGraph();
  grRelWF->SetName("grRelWF");
  grRelWF->Set(NumrBins);
  TGraph *grImWF= new TGraph();
  grImWF->SetName("grImWF");
  grImWF->Set(NumrBins);
  TGraph *grDensityWF= new TGraph();
  grDensityWF->SetName("grDensityWF");
  grDensityWF->Set(NumrBins);
  TLegend *leg1= new TLegend(0.6, 0.4, 0.85, 0.55);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.040);
  leg1->SetLineColor(0);
  leg1->AddEntry(grRelWF, "#Rgothic [#Phi]", "pef");
  TLegend *leg2= new TLegend(0.6, 0.3, 0.85, 0.5);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.040);
  leg2->SetLineColor(0);
  leg2->AddEntry(grImWF, "#Jgothic [#Phi]", "pef");
  TLegend *leg3= new TLegend(0.6, 0.3, 0.85, 0.5);
  leg3->SetFillStyle(0);
  leg3->SetTextSize(0.040);
  leg3->SetLineColor(0);
  leg3->AddEntry(grRelWF, "#Rgothic [#Phi]", "pef");
  leg3->AddEntry(grImWF, "#Jgothic [#Phi]", "pef");
  TLegend *leg4= new TLegend(0.6, 0.5, 0.85, 0.7);
  leg4->SetFillStyle(0);
  leg4->SetTextSize(0.040);
  leg4->SetLineColor(0);
  leg4->AddEntry(grDensityWF, "#Phi^{2} Prob density", "pef");
  // leg2->AddEntry(grRelWF, "Chiral EFT", "pef");
  for (unsigned uMom= 0; uMom < NumrBins; uMom++) {
    grRelWF->SetPoint(uMom, uMom, Evaluate_d_u(uMom));  //
    grImWF->SetPoint(uMom, uMom, Evaluate_d_w(uMom));  //
    grDensityWF->SetPoint(
        uMom, uMom,
        pow(Evaluate_d_u(uMom), 2.0) + pow(Evaluate_d_w(uMom), 2.0));
    std::cout << M1->GetMomentum(uMom) << "    " << Evaluate_d_u(uMom) << "    "
        << Evaluate_d_w(uMom) << std::endl;

  }

  grDensityWF->SetTitle("; #it{r} (fm); ");
  // grDensityWF->GetYaxis()->SetRangeUser(-0.01, 0.35);
 // grDensityWF->SetTitle("k = 20MeV, #Lambda = 800");
//  grImWF->SetTitle("k = 20MeV, #Lambda = 800");
 // grRelWF->SetTitle("k = 20MeV, #Lambda = 800");

  grDensityWF->SetMarkerColor(kRed + 2);
  grDensityWF->SetLineColor(kRed + 2);
  grDensityWF->GetXaxis()->SetTitleSize(0.05);
  //grDensityWF->GetYaxis()->SetTitleSize(0.05);
  //grDensityWF->GetXaxis()->SetLabelSize(0.05);
  grDensityWF->GetXaxis()->SetTitleOffset(0.6);
  grDensityWF->GetYaxis()->SetLabelSize(0.03);
  //grDensityWF->GetYaxis()->SetTitleOffset(0.6);

  grImWF->SetMarkerColor(65 + 2);
  grImWF->SetLineColor(65 + 2);
  grImWF->SetTitle("; #it{r} (fm); ");
  grImWF->GetXaxis()->SetTitleSize(0.05);
  grImWF->GetYaxis()->SetTitleSize(0.05);
  grImWF->GetXaxis()->SetTitleOffset(0.6);
  grImWF->GetYaxis()->SetTitleOffset(0.6);
  grRelWF->SetMarkerColor(kBlue + 2);
  grRelWF->SetLineColor(kBlue + 2);
  grRelWF->SetTitle("; #it{r} (fm); ");
  grRelWF->GetXaxis()->SetTitleSize(0.05);
  grRelWF->GetYaxis()->SetTitleSize(0.05);
  grRelWF->GetXaxis()->SetTitleOffset(0.6);
  grRelWF->GetYaxis()->SetTitleOffset(0.6);
  /*
  p1->cd();
  grRelWF->Draw("ALP");
  leg1->Draw("same");

  p2->cd();
  grImWF->Draw("ALP");
  leg2->Draw("same");
  p3->cd();
  grImWF->Draw("ALP");
  grRelWF->Draw("same");
  leg3->Draw("same");
  p4->cd();
  grDensityWF->Draw("ALP");
  leg4->Draw("same");*/

 //grDensityWF->Draw("ALP");
  WFplot->cd(1);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.1);
  gPad->SetRightMargin(0.14);
  gPad->SetLeftMargin(0.12);
  grRelWF->Draw("ALP");
  leg1->Draw("same");
  WFplot->cd(2);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.1);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.1);

  grImWF->Draw("ALP");
  leg2->Draw("same");
  WFplot->cd(3);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.14);
  gPad->SetLeftMargin(0.12);
  grImWF->Draw("ALP");
  grRelWF->Draw("same");
  leg3->Draw("same");
  WFplot->cd(4);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.1);
  grDensityWF->Draw("ALP");
  leg4->Draw("same");
  WFplot->SaveAs(nameWFplot);
  WFplot->Close();
 // c1->SaveAs(nameWFplot);
 // c1->Close();
  delete M1;
  // delete M2;
  // delete M3;

}
