#include "Model.h"
Model::Model()
    : frmin(0.0),
      frmax(0.0),
      fWFGauss(0),
      fWFTwoGauss(0),
      fWFHulthen(0),
      MomBin(nullptr),
      MomBinCenter(nullptr),
      // MomBinCenter(nullptr),
      rWF(nullptr),
      kWFErr(nullptr),
      fWFchiEFT(0) {
  NumMomBins= 0;
}
Model::~Model() {
  if (MomBin) {
    delete[] MomBin;
    MomBin= NULL;
  }
  if (MomBinCenter) {
    delete[] MomBinCenter;
    MomBinCenter= NULL;
  }
}

void Model::DelAllMom() {
  DelMom();
}
double Model::GetMomentum(const unsigned &WhichMomBin) const {
  if (NumMomBins <= WhichMomBin)
    return 0;
  return MomBinCenter[WhichMomBin];
}

double Model::GetMomBinLowEdge(const unsigned &WhichMomBin) const {
  if (NumMomBins < WhichMomBin)
    return 0;
  return MomBin[WhichMomBin];
}

double Model::GetMomBinUpEdge(const unsigned &WhichMomBin) const {
  if (NumMomBins <= WhichMomBin)
    return 0;
  return MomBin[WhichMomBin + 1];
}

void Model::SetMomBins(const unsigned &nummombins, const double *mombins,
                       const double *bincenter) {
  if (!nummombins) {
    if (Notifications >= nError)
      printf(
          "\033[1;31mERROR:\033[0m Bad input in CATS::SetMomBins(const unsigned& nummombins, const double* mombins)\n");
    return;
  }
  if (!mombins) {
    if (Notifications >= nError)
      printf(
          "\033[1;31mERROR:\033[0m Bad input in CATS::SetMomBins(const unsigned& nummombins, const double* mombins)\n");
    return;
  }
  //check if the momentum bins set are the same as before. If yes, change nothing
  if (nummombins == NumMomBins) {
    bool SameBinning= true;
    for (unsigned uMomBin= 0; uMomBin < NumMomBins; uMomBin++) {
      SameBinning*= (mombins[uMomBin] == MomBin[uMomBin]);
    }
    if (SameBinning)
      return;
  }
  if (nummombins != NumMomBins || !MomBin) {
    if (MomBin) {
      delete[] MomBin;
      MomBin= NULL;
    }
    if (MomBinCenter) {
      delete[] MomBinCenter;
      MomBinCenter= NULL;
    }
    MomBin= new double[nummombins + 1];
    MomBinCenter= new double[nummombins];
    DelAllMom();
    NumMomBins= nummombins;
  }

  for (unsigned uBin= 0; uBin <= NumMomBins; uBin++) {
    MomBin[uBin]= mombins[uBin];
    if (uBin) {
      if (bincenter)
        MomBinCenter[uBin - 1]= bincenter[uBin - 1];
      else
        MomBinCenter[uBin - 1]= 0.5 * (mombins[uBin - 1] + mombins[uBin]);
    }

    if (MomBin[uBin] < 0) {
      if (Notifications >= nError) {
        printf(
            "\033[1;31mERROR:\033[0m CATS::SetMomBins(const unsigned& nummombins, const double* mombins)\n");
        printf("         The momentum should be positive!\n");
      }
      return;
    }
  }
}
void Model::SetMomBins(const unsigned &nummombins, const double &MinMom,
                       const double &MaxMom) {
  if (!nummombins) {
    if (Notifications >= nError)
      printf(
          "\033[1;31mERROR:\033[0m Bad input in CATS::SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom)\n");
    return;
  }
  if (MinMom > MaxMom) {
    if (Notifications >= nError)
      printf(
          "\033[1;31mERROR:\033[0m Bad input in CATS::SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom)\n");
    return;
  }
  if (MinMom == MaxMom && nummombins != 1) {
    if (Notifications >= nError)
      printf(
          "\033[1;31mERROR:\033[0m Bad input in CATS::SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom)\n");
    return;
  }
  //check if the momentum bins set are the same as before. If yes, change nothing
  double BinWidth= (MaxMom - MinMom) / double(nummombins);
  if (nummombins == NumMomBins && MomBin[0] == MinMom
      && MomBin[nummombins] == MinMom + double(nummombins) * BinWidth)
    return;
  if (nummombins != NumMomBins || !MomBin) {
    if (MomBin) {
      delete[] MomBin;
      MomBin= NULL;
    }
    if (MomBinCenter) {
      delete[] MomBinCenter;
      MomBinCenter= NULL;
    }
    MomBin= new double[nummombins + 1];
    MomBinCenter= new double[nummombins];
    DelAllMom();
    NumMomBins= nummombins;
  }

  for (unsigned uBin= 0; uBin <= NumMomBins; uBin++) {
    MomBin[uBin]= MinMom + double(uBin) * BinWidth;
    if (uBin != NumMomBins)
      MomBinCenter[uBin]= MinMom + double(uBin) * BinWidth + 0.5 * BinWidth;
  }
}
void Model::DelMom() {
  if (rWF) {
    delete[] rWF;
    rWF= NULL;
    delete[] kWFErr;
    kWFErr= NULL;
  }
}

double Model::GetWaveFunction(const unsigned &WhichMomBin) {
  if (WhichMomBin >= NumMomBins) {
    std::cout << "return is shitty\n";
    return 0;
  }
  if (fWFGauss) {
    return Wavefun.WFGauss(MomBinCenter[WhichMomBin]);
  }

  if (fWFTwoGauss) {

  }
  if (fWFHulthen) {
    return Wavefun.WFHulthen(MomBinCenter[WhichMomBin]);
  }
  if (fWFchiEFT) {
    return Wavefun.WaveChiEFT(MomBinCenter[WhichMomBin]);
  } else {
    std::cout << " No Wavefunction to draw r\n ";
  }
}
//This gives you B2 using Gaussian wave function
double Model::GetB2_Classical(const unsigned &WhichMomBin) {
  double Rval;
  Rval= MomBinCenter[WhichMomBin];
  double Norm= 7.9775;//at R=6 N= 7.97747772984 and at R = 10,N= 8.4343411;  //this normalization is always taken with respect to B2 quantum  at large R
  double Unit_con_Factor= 1 / 5.068 * 1 / 5.068 * 1 / 5.068;  //fm^{-1} to GeV
  return Norm * 1 / pow(Rval, 3.0) * Unit_con_Factor;
}  //This gives you B2 using Gaussian wave function
double Model::GetB2_kfir(const unsigned &WhichMomBin) {
  double Rval;
  Rval= MomBinCenter[WhichMomBin];
  double par1= 0.938;  //GeV
  double par2= 3.2;  //fm
  double Unit_con_Factor= 1 / 5.068 * 1 / 5.068 * 1 / 5.068;  //fm^{-1} to GeV
  return (3 * pow(TMath::Pi(), 3. / 2.))
      / (2. * par1 * pow((Rval * Rval + par2 * par2 / 4.), 3. / 2.))
      * Unit_con_Factor;
}

double Model::GetB2_HulthenKfir(const unsigned &WhichMomBin) {
  double Rval;
  Rval= MomBinCenter[WhichMomBin];
  float a= 0.2;  // alpha
  float b= 1.56;  // Beta
  float N= 3.0 / 0.94;
  double Unit_con_Factor= 1 / 5.068 * 1 / 5.068 * 1 / 5.068;  //fm^{-1} to GeV

  return N * TMath::Pi() * TMath::Pi() / pow(Rval, 2.0) * (a * b * (a + b))
      / ((a - b) * (a - b))
      * (TMath::Exp(4 * a * a * Rval * Rval) * TMath::Erfc(2 * a * Rval)
          - 2 * TMath::Exp((a + b) * (a + b) * Rval * Rval)
              * TMath::Erfc((a + b) * Rval)
          + TMath::Exp(4 * b * b * Rval * Rval) * TMath::Erfc(2 * b * Rval))*Unit_con_Factor;

}
//This gives you dB2(r,q,R) using Hulthen wave function
double Model::dB2Hulthen(double *x, double *par) {
  double qval= x[0];
  double Rval;
  Rval= par[0];
  double r= x[1];
  float a= 0.2;  // alpha
  float b= 1.56;  // Beta
  float N= 1.5 / 0.94;
  double Unit_con_Factor= 1 / 5.068 * 1 / 5.068 * 1 / 5.068;  //fm^{-1} to GeV
  return N * 16 * TMath::Pi() * TMath::Pi() * qval * qval
      * TMath::Exp(-pow(Rval * qval, 2.0)) * TMath::Sin(qval * r) / (qval * r)
      * (a * b * (a + b) / (2 * TMath::Pi() * (a - b) * (a - b)))
      * (TMath::Exp(-a * r) - TMath::Exp(-b * r))
      * (TMath::Exp(-a * r) - TMath::Exp(-b * r)) * Unit_con_Factor;
}

//Hulthen wave function does not converge at very low R and therefore the integration needs to be done numerically
//I did used cheap TF2 function to have numerical integration
//the uncertainty due to this trick is <<<< than the existing uncertainty in Coalescence Model
//Integral r limits are chosen such that the error in B2 is 10^{-6} order which is very good approximation
//considering that fact that other source of error are more than 10%

double Model::GetB2Hulthen(const unsigned &WhichMomBin) {
  double B2_HulthenValue;
  double *x;
  double *par;
  float q_limit= 10.0;  // it's in fm^{-1}, We choose this limit after checking the D[q] vs q ...the error you introduce is 10^{-6} order.
  TF2 *B2_Hulthen= new TF2("B2_Hulthen1", dB2Hulthen, 0.0001, 100.0, 0.00,
                           10000000, 1);
  B2_Hulthen->SetParameter(0, MomBinCenter[WhichMomBin]);
  B2_HulthenValue= B2_Hulthen->Integral(0.01, q_limit, 0.0, 324.0);
  return B2_HulthenValue;
}

double Model::dB2TwoGaussian(double *x, double *par) {
  double qval= x[0];
  double Rval;
  Rval= par[0];
  double r= x[1];
  double del= 0.581;  //del
  double d1= 3.979;  //d1
  double d2= 0.890;  //d2
  float N= 1.5 / 0.94;
  double Unit_con_Factor= 1 / 5.068 * 1 / 5.068 * 1 / 5.068;  //fm^{-1} to GeV
  return N * 16 * pow(TMath::Pi(), 0.5) * qval * r
      * TMath::Exp(-pow(Rval * qval, 2.0)) * TMath::Sin(qval * r)
      * (del / pow(d1, 3.0) * TMath::Exp(-pow(r / d1, 2.0))
          + (1 - del) / d2 * TMath::Exp(-pow(r / d2, 2.0))) * Unit_con_Factor;
}
double Model::GetB2TwoGaussian(const unsigned &WhichMomBin) {
  double B2_TwoGaussianValue;
  double *x;
  double *par;
  float q_limit= 10.0;  // it's in fm^{-1}, We choose this limit after checking the D[q] vs q ...the error you introduce is 10^{-6} order.
  TF2 *B2TwoGaussian= new TF2("dB2TwoGaussian", dB2TwoGaussian, 0.0000, 100.0,
                              0.00, 10000000, 1);
  B2TwoGaussian->SetParameter(0, MomBinCenter[WhichMomBin]);
  B2_TwoGaussianValue= B2TwoGaussian->Integral(0.00, q_limit, 0.0, 324.0);
  return B2_TwoGaussianValue;
}

double Model::GetEvaluate_d_u(double r) {
  return Wavefun.uWaveChiEFT(r);
}
;
double Model::GetEvaluate_d_w(double r) {
  return Wavefun.wWaveChiEFT(r);
}
void Model::LoadWF() {
  Wavefun.LoadTheWF();
}
double Model::dB2ChiEFT(double *x, double *par) {
  float qval= x[0];
  double R= par[0];
  float r= x[1];
  float N= 1.5 / 0.94;
  double Unit_con_Factor= 1 / 5.068 * 1 / 5.068 * 1 / 5.068;  //fm^{-1} to GeV
  Model *myodel= new Model();
  double u= myodel->GetEvaluate_d_u(r);
  double w= myodel->GetEvaluate_d_w(r);
  return N * 4 * TMath::Pi() * qval * TMath::Exp(-pow((R * qval), 2.0))
      * (u * u + w * w) * TMath::Sin(qval * r) / r * Unit_con_Factor;
}
double Model::GetB2ChiEFT(const unsigned &WhichMomBin) {
  double B2_ChiEFTValue;
  double *x;
  double *par;
  float q_limit= 3.07;  //this limit here is because of ChiEFT cutt off
  TF2 *B2_ChiEFT1= new TF2("B2_ChiEFT1", dB2ChiEFT, 0.0, 50.5, 0.05, 14.5, 1);
  B2_ChiEFT1->SetParameter(0, MomBinCenter[WhichMomBin]);
  B2_ChiEFTValue= B2_ChiEFT1->Integral(0.0, q_limit, 0.01, 14.05);
  return B2_ChiEFTValue;
}

